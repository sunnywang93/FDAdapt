#' Estimate global adaptive bandwidth for covariance function
#'
#' `estimate_bandwidth_covariance_norm` estimates the global adaptive bandwidth
#' in terms of integrated squared risk, as opposed to
#' `estimate_bandwidth_covariance`, which estimates a pointwise bandwidth.
#'
#' @param curves List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param grid_bandwidth Grid of points to estimate bandwidth.
#' @param grid_smooth Grid of points to smooth curves.
#' @param k0 Minimum number of points curves must possess to be used in
#' estimation.
#' @param grid_param Grid of points to estimate parameters.
#' @param sigma Noise level, if known. Defaults to NULL, which estimates it.
#' @param mu0 Density lower bound for time points. Defaults to NULL, which
#' estimates it.
#' @returns A list with two elements, containing
#' * A global bandwidth.
#' * A list, containing
#'   - **$moments** Estimated moments.
#'   - **$constants** Estimated Hölder constants `H(t)` and `L(t)`.
#'   - **$kernel_int** Integrated constant of kernel.
#'   - **$N_gamma_ts** Harmonic mean of the number of points used for
#'   smoothing each curve at `(t|s)`.
#'   - **$N_gamma_st** Harmonic mean of the number of points used for
#'   smoothing each curve at `(s|t)`.
#'   - **$WN** Total number of curves selected for estimation.
#' @export

estimate_bandwidth_covariance_norm <- function(curves, grid_bandwidth,
                                               grid_smooth, k0,
                                               grid_tibble)
{

  cst_kernel <- 1.5 * (1 / (1 + 2 * grid_tibble$H) - 1 / (3 + 2 * grid_tibble$H))

  M_length <- curves |> purrr::map_dbl(~(length(.x$t)))
  M_length_cum <- cumsum(M_length)

  vec <- seq(1, max(M_length_cum))
  vec_cut <- cut(vec, breaks = c(1, M_length_cum), labels = FALSE,
                 include.lowest = TRUE)
  #seq of indices
  Mi_vec <- split(vec, vec_cut) |> unname()

  T_diff <- lapply(curves, function(x) x$t) |>
    unlist() |>
    outer(grid_smooth, FUN = "-")

  indic <- T_diff |>
    abs() |>
    outer(grid_bandwidth, FUN = "<=") * 1

  #dim G x N x H
  ws <- sapply(Mi_vec, function(Mi) {
    indic[Mi,,] |> colSums(na.rm = TRUE) |>
      (\(x) x >= k0)() * 1
  }, simplify = "array") |>
    aperm(c(1, 3, 2))

  #dim G x G x N x H
  wst <- apply(ws, MARGIN = c(2, 3), function(x) x %*% t(x)) |>
    array(dim = c(length(grid_smooth),
                  length(grid_smooth),
                  length(curves),
                  length(grid_bandwidth)))


  #dim G x G x H
  WN <- wst |> aperm(c(3, 4, 1, 2)) |>
    colSums(na.rm = TRUE) |>
    aperm(c(2, 3, 1))

  Wm_t <- outer(T_diff, grid_bandwidth, FUN = "/") |>
    epa_kernel()
  #N x G x H
  max_Wm_t <- sapply(Mi_vec, function(Mi) {
    Wm_t[Mi,,] |> apply(c(2,3), FUN = max, na.rm = TRUE)
  }, simplify = "array")

  sum_Wm_t <- sapply(Mi_vec, function(Mi) {
    Wm_t[Mi,,] |> colSums(na.rm = TRUE)
  }, simplify = "array")

  max_Wm_t_norm <- (max_Wm_t / sum_Wm_t)

  rm(max_Wm_t, sum_Wm_t, Wm_t)

  max_Wm_t_norm <- aperm(max_Wm_t_norm, c(1, 3, 2))

  Ngamma_ts <- array(max_Wm_t_norm,
                     dim = c(dim(max_Wm_t_norm), length(grid_smooth))) |>
    (\(x) wst / aperm(x, c(4, 1, 2, 3)))() |>
    (\(Ni) aperm(wst / Ni, c(3,1,2,4)))() |>
    colSums(na.rm = TRUE) |>
    (\(x) 1 / (x / WN**2))()

  Ngamma_st <- aperm(Ngamma_ts, c(2, 1, 3))

  var_ <- estimate_variance_curves(data = curves, params = grid_tibble,
                                   grid_smooth = grid_smooth)

  q1_ts <- outer(cst_kernel * grid_tibble$L**2, 2 * var_$EXt2)

  h_alpha_t <- t(outer(grid_bandwidth, 2 * grid_tibble$H, FUN = "^"))

  ones <- rep(1, length(grid_smooth))

  q1_ts_term <- sapply(seq_along(grid_bandwidth), function(h)
    outer(h_alpha_t[, h], ones) * q1_ts, simplify = "array")

  q1_st_term <- aperm(q1_ts_term, c(2, 1, 3))


  qq2 <- matrix(var_$EXt2, ncol = length(grid_smooth),
                nrow = length(grid_smooth))

  #q2(t|s) = q2(s|t) = q2
  q2 <- max(grid_tibble$sigma**2) * (qq2 + t(qq2))

  q2_ts_term <- array(q2, dim = c(length(grid_smooth),
                                  length(grid_smooth),
                                  length(grid_bandwidth))) /
    Ngamma_ts

  q2_st_term <- array(q2, dim = c(length(grid_smooth),
                                  length(grid_smooth),
                                  length(grid_bandwidth))) /
    Ngamma_st

  q3 <- var_$varXtXs / 2

  q3_term <- replicate(length(grid_bandwidth), q3) * (1/WN - 1/length(curves))

  #risk(t|s) is the transpose of risk(s|t)
  #risk <- risk_s + aperm(risk_s, c(2, 1, 3)) + 2 * q3_term
  risk <- q1_ts_term + q1_st_term +
    q2_ts_term + q2_st_term +
    2 * q3_term

  integrated_risk <- apply(risk, c(2, 3),
                           function(t) pracma::trapz(grid_smooth, t)) |>
    apply(2, function(s) pracma::trapz(grid_smooth, s))

  min_h_index <- which.min(integrated_risk)

  h_opt <- grid_bandwidth[min_h_index]

  M_avg <- purrr::map_dbl(curves, ~length(.x$t)) |> mean()

  N_effective <- mean(WN[,, min_h_index]) * M_avg

  h_constant <- log(N_effective)**(abs(log(h_opt) / log(N_effective)))

  list(bandwidth = h_constant * h_opt,
       params = list(moments = var_,
                     constants = grid_tibble,
                     kernel_int = cst_kernel,
                     Ngamma_ts = Ngamma_ts,
                     Ngamma_st = Ngamma_st,
                     WN = WN))

}

#' Smooth curves using adaptive global bandwidths for covariance function
#' estimation (in terms of `L2` norm).
#'
#' `smooth_curves_covariance_norm` smooths each curve using the bandwidth
#' estimated from `estimate_bandwidth_covariance_norm`, for the purpose
#' of covariance function estimation.
#'
#' @param curves List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param grid_bandwidth Grid of points to estimate bandwidth.
#' @param grid_smooth Grid of points to smooth curves.
#' @param k0 Minimum number of points curves must possess to be used in
#' estimation.
#' @param grid_param Grid of points to estimate parameters.
#' @param sigma Noise level, if known. Defaults to NULL, which estimates it.
#' @param mu0 Density lower bound for time points. Defaults to NULL, which
#' estimates it.
#' @returns A list with two elements, containing
#'   - **$smoothed_product** List containing the smoothed product `XtXs`.
#'   - **$smoothed_curves** List containing the smoothed curves `Xt`.
#'   - **$bandwidth** Global bandwidth used to smooth each curve.
#'   - **$params** A list, containing
#'     - **$moments** Estimated moments.
#'     - **$constants** Estimated Hölder constants `H(t)` and `L(t)`.
#'     - **$kernel-int** Integrated constant of kernel.
#'     - **$N_gamma_ts** Harmonic mean of the number of points used for
#'   smoothing each curve at `(t|s)`.
#'     - **$N_gamma_st** Harmonic mean of the number of points used for
#'   smoothing each curve at `(s|t)`.
#'     - **$WN** Total number of curves selected for estimation.
#' @examples
#' smooth_curves_covariance(curves = curves_list,
#' grid_bandwidth = seq(0, 1, length.out = 151),
#' grid_smooth = seq(0, 1, length.out = 101),
#' k0 = 1,
#' grid_param = seq(0, 1, length.out = 20),
#' sigma = 0.1, mu0 = 1, nvalues = 10)
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#' estimation of irregular mean and covariance functions.
#' @export


smooth_curves_covariance_norm <- function(curves, grid_bandwidth,
                                          grid_smooth, k0, grid_tibble)
{

  bandwidth_list <- estimate_bandwidth_covariance_norm(curves,
                                                       grid_bandwidth,
                                                       grid_smooth,
                                                       k0,
                                                       grid_tibble)

  M_length <- curves |> purrr::map_dbl(~(length(.x$t)))

  #dim Mi x G
  Wm_num <- purrr::map(curves,
                       ~epa_kernel(abs(outer(.x$t, grid_smooth, FUN = "-")) /
                                     bandwidth_list$bandwidth))

  Wm_denom <- purrr::map(Wm_num, ~colSums(.x, na.rm = TRUE))

  Wm <- purrr::map2(Wm_num, Wm_denom,
                    ~sweep(.x, MARGIN = 2, STATS = .y, FUN = "/") |>
                      (\(x) replace(x, is.nan(x), 0))())

  Xt <- purrr::map2(Wm, curves, ~t(.x)%*%.y$x)

  XtXs <- purrr::map(Xt, ~.x %*% t(.x))

  list(smoothed_product = XtXs,
       smoothed_curves = Xt,
       bandwidth = bandwidth_list$bandwidth,
       params = bandwidth_list$params)
}

#' Estimate covariance function using adaptive global bandwidth
#'
#' `covariance_norm` estimates the covariance function using curves that
#' have been smoothed using an adaptive global bandwidth.
#'
#' @param curves List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param grid_bandwidth Grid of points to estimate bandwidth.
#' @param grid_smooth Grid of points to smooth curves.
#' @param k0 Minimum number of points curves must possess to be used in
#' estimation.
#' @param grid_param Grid of points to estimate parameters.
#' @param sigma Noise level, if known. Defaults to NULL, which estimates it.
#' @param mu0 Density lower bound for time points. Defaults to NULL, which
#' estimates it.
#' @returns A list containing
#' - **$cov** Estimated covariance function.
#' - **$bandwidth** Global bandwidth used to smooth each curve.
#' - **$constants** Tibble containing Hölder constants.
#' - **$kernel_int** Integrated kernel.
#' - **$Ngamma_ts** Harmonic mean of the number of points used to smooth
#' each curve at `t|s`.
#' - **$Ngamma_st** Transpose of `Ngamma_ts`.
#' - **$WN** Total number of points selected for estimation.
#' @examples
#' covariance_ll(curves = curves_list,
#' grid_bandwidth = seq(0, 1, length.out = 151),
#' grid_smooth = seq(0, 1, length.out = 101),
#' k0 = 1,
#' grid_param = seq(0, 1, length.out = 20),
#' sigma = 0.1, mu0 = 1, nvalues = 10)
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#' estimation of irregular mean and covariance functions.
#' @export


covariance_norm <- function(curves, grid_bandwidth, grid_smooth,
                            k0, params) {

  m <- purrr::map_dbl(curves, ~length(.x$t)) |> mean()

  interp_grid <- c(0, params$t, 1)

  var_interp <- as.matrix(params[, c("H", "L", "sigma", "mu0")])

  interp_mat <- apply(var_interp, MARGIN = 2,
                      function(x) pracma::interp1(x = interp_grid,
                                                  y = c(x[1], x, x[nrow(var_interp)]),
                                                  xi = grid_smooth))

  grid_tibble <- tibble::as_tibble(cbind(t = grid_smooth, interp_mat))

  prod_ <- smooth_curves_covariance_norm(curves, grid_bandwidth,
                                         grid_smooth, k0, grid_tibble)

  indic <- purrr::map(curves,
                   ~(abs(outer(.x$t, grid_smooth, "-")) <= prod_$bandwidth) * 1)

  wt <- purrr::map(indic, ~(colSums(.x, na.rm = TRUE) >= k0) * 1)

  WN_t <- Reduce('+', wt)

  wtws <- purrr::map(wt, ~.x %*% t(.x))

  WN <- Reduce('+', wtws)

  prod_sum <- purrr::map2(wtws, prod_$smoothed_product, ~(.x * .y)) |>
    (\(x) Reduce('+', x))()

  gamma <- (1/WN) * prod_sum

  mu_t <- purrr::map2(wt, prod_$smoothed_curves, ~.x * .y) |>
    (\(x) Reduce('+', x))() / WN_t

  mu <- mu_t %*% t(mu_t)

  Gamma <- gamma - mu

  #final step: replace diagonal band
  for (t in 1:ncol(Gamma)) {
    s <- 1
    current_cov <- Gamma[s, t - s + 1]
    while (s <= (t - s + 1)) {
      if (abs(grid_smooth[s] - grid_smooth[t - s + 1]) >
          prod_$bandwidth) {
        current_cov <- Gamma[s, t - s + 1]
      } else {
        Gamma[s, t - s + 1] <- current_cov
      }
      s <- s + 1
    }
  }
  #loop over lower anti-diagonal
  for (s in 1:nrow(Gamma)) {
    t <- ncol(Gamma)
    current_cov <- Gamma[ncol(Gamma) + s - t, t]
    while (t >= (ncol(Gamma) + s - t)) {
      if (abs(grid_smooth[ncol(Gamma) + s - t] - grid_smooth[t]) >
          prod_$bandwidth) {
        current_cov <- Gamma[ncol(Gamma) + s - t, t]
      } else {
        Gamma[ncol(Gamma) + s - t, t] <- current_cov
      }
      t <- t - 1
    }
  }
  list(cov = Gamma,
       bandwidth = prod_$bandwidth,
       constants = prod_$params$constants,
       kernel_int = prod_$params$kernel_int,
       Ngamma_ts = prod_$params$Ngamma_ts,
       Ngamma_st = prod_$params$Ngamma_st,
       variance = prod_$params$moments$varXt,
       variance_prod = prod_$params$moments$varXtXs,
       moment2 = prod_$params$moments$EXt2,
       moment2_prod = prod_$params$moments$EXtXs2,
       WN = prod_$params$WN)
}



