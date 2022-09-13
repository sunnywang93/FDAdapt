#' Estimate adaptive bandwidths for covariance function
#'
#' `estimate_bandwidth_covariance` estimates the bandwidth at every point
#' `(t,s)` that is used to smooth each curve for the purpose of
#' covariance function estimation.
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
#' * A matrix of bandwidth for every point `(t,s)`.
#' * A list, containing
#'   - **$moments** Estimated moments.
#'   - **$constants** Estimated Hölder constants `H(t)` and `L(t)`.
#'   - **$kernel-int** Integrated constant of kernel.
#'   - **$N_gamma_ts** Harmonic mean of the number of points used for
#'   smoothing each curve at `(t|s)`.
#'   - **$N_gamma_st** Harmonic mean of the number of points used for
#'   smoothing each curve at `(s|t)`.
#'   - **$WN** Total number of curves selected for estimation.
#' @examples
#' estimate_bandwidth_covariance(curves = curves_list,
#' grid_bandwidth = seq(0, 1, length.out = 151),
#' grid_smooth = seq(0, 1, length.out = 101),
#' k0 = 1,
#' grid_param = seq(0, 1, length.out = 20),
#' sigma = 0.1, mu0 = 1, nvalues = 10)
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#' estimation of irregular mean and covariance functions.
#' @export


estimate_bandwidth_covariance <- function(curves, grid_bandwidth,
                          grid_smooth, k0, grid_param, sigma = NULL,
                          mu0 = NULL)
{

  params <- estimate_holder_const(curves, sigma = sigma, mu0 = mu0,
                                  grid_estim = grid_param)

  H_grid_smooth <- pracma::interp1(x = c(0, grid_param, 1),
                                   y = c(params$H[1], params$H,
                                         params$H[length(params$H)]),
                                   xi = grid_smooth, method = "linear")

  L_grid_smooth <- pracma::interp1(x = c(0, grid_param, 1),
                                   y = c(params$L[1], params$L,
                                         params$L[length(params$L)]),
                                   xi = grid_smooth, method = "linear")


  grid_tibble <- tibble::tibble(t = grid_smooth,
                                H = H_grid_smooth,
                                L = L_grid_smooth,
                                sigma = max(sigma))

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
  max_Wm_t_norm[is.nan(max_Wm_t_norm)] <- 0

  max_Wm_t_norm <- aperm(max_Wm_t_norm, c(1, 3, 2))

  gamma_sum_ts <- sapply(seq_along(grid_bandwidth), function(h) {
    max_Wm_t_norm[,,h] %*% t(ws[,,h])
  }, simplify = "array")

  #t are columns, s are rows, $N_{\Gamma}(t|s)$
  Ngamma_ts <- 1 / (1/WN^2 * gamma_sum_ts)
  Ngamma_ts[!is.finite(Ngamma_ts)] <- 10^(-100)

  Ngamma_st <- aperm(Ngamma_ts, c(2, 1, 3))
  Ngamma_st[!is.finite(Ngamma_st)] <- 10^(-100)


  var_ <- estimate_variance_curves(data = curves, params = params,
                                   grid_smooth = grid_smooth,
                                   sigma = grid_tibble$sigma,
                                   mu0 = unique(params$mu0))


  q1_ts <- sapply(seq_along(grid_smooth), function(s) {
    sqrt(2 * var_$EXt2[s] * cst_kernel / factorial(floor(grid_tibble$H))^2 *
           grid_tibble$L)
  })

  q1_st <- t(q1_ts)

  ones <- rep(1, length(grid_smooth))

  h_alpha_t <- sapply(grid_bandwidth, function(h) {
    h^(2*grid_tibble$H)
  })

  q1_ts_term <- sapply(seq_along(grid_bandwidth), function(h)
    outer(h_alpha_t[, h], ones) * q1_ts, simplify = "array")

  q1_st_term <- aperm(q1_ts_term, c(2, 1, 3))

  qq2 <- matrix(var_$EXt2, ncol = length(grid_smooth),
                nrow = length(grid_smooth))

  #q2(t|s) = q2(s|t) = q2
  q2 <- max(sigma^2) * (qq2 + t(qq2))

  q2_ts_term <- array(q2, dim = c(length(grid_smooth),
                               length(grid_smooth),
                               length(grid_bandwidth))) /
    Ngamma_ts

  q2_st_term <- array(q2, dim = c(length(grid_smooth),
                                  length(grid_smooth),
                                  length(grid_bandwidth))) /
    Ngamma_st

  q3 <- sqrt(var_$varXtXs / 2)

  q3_term <- replicate(length(grid_bandwidth), q3^2) * (1/WN - 1/length(curves))

  #risk(t|s) is the transpose of risk(s|t)
  #risk <- risk_s + aperm(risk_s, c(2, 1, 3)) + 2 * q3_term
  risk <- q1_ts_term + q1_st_term +
    q2_ts_term + q2_st_term +
    2 * q3_term

  min_h_index <- apply(risk, MARGIN = c(1,2), which.min)

  list(bw_matrix = apply(min_h_index, 2, function(id) grid_bandwidth[id]),
       params = list(moments = var_,
                     constants = grid_tibble,
                     kernel_int = cst_kernel,
                     Ngamma_ts = Ngamma_ts,
                     Ngamma_st = Ngamma_st,
                     WN = WN))

}

#' Smooth curves using adaptive bandwidths for covariance function
#' estimation
#'
#' `smooth_curves_covariance` smooths each curve using the bandwidth
#' estimated from `estimate_bandwidth_covariance`, for the purpose
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
#' * The product of curves `XtXs`.
#' * A list, with
#'   - A matrix of bandwidth for every point `(t,s)`.
#'   - A list, containing
#'   - **$moments** Estimated moments.
#'   - **$constants** Estimated Hölder constants `H(t)` and `L(t)`.
#'   - **$kernel-int** Integrated constant of kernel.
#'   - **$N_gamma_ts** Harmonic mean of the number of points used for
#'   smoothing each curve at `(t|s)`.
#'   - **$N_gamma_st** Harmonic mean of the number of points used for
#'   smoothing each curve at `(s|t)`.
#'   - **$WN** Total number of curves selected for estimation.
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


smooth_curves_covariance <- function(curves, grid_bandwidth,
                                     grid_smooth, k0, grid_param,
                                     sigma, mu0)
  {

  bandwidth_list <- estimate_bandwidth_covariance(curves,
                                                  grid_bandwidth,
                                                  grid_smooth,
                                                  k0,
                                                  grid_param,
                                                  sigma,
                                                  mu0)

  bandwidth_matrix <- bandwidth_list$bw_matrix

  M_length <- curves |> purrr::map_dbl(~(length(.x$t)))

  #dim Mi x G
  Wm_num <- purrr::map(curves, ~outer(.x$t, grid_smooth, "-")) |>
    lapply(function(Ti) {
      sapply(seq_along(grid_smooth), function(t) {
        outer(Ti[, t], bandwidth_matrix[, t], FUN = "/")
      }, simplify = "array") |> epa_kernel()
    })

  Wm_denom <- purrr::map(Wm_num, ~colSums(.x, na.rm = TRUE)) |>
    purrr::map2(M_length, ~replicate(.y, .x)) |>
    purrr::map(~aperm(.x, c(3, 1, 2)))

  Wm <- purrr::map2(Wm_num, Wm_denom, ~(.x / .y) |>
                      (\(x) replace(x, is.nan(x), 0))())

  XtXs <- purrr::map2(curves, Wm, function(i, w) {
    sapply(seq_along(grid_smooth), function(s) {
      t(w[,,s]) %*% i$x
    })
  }) |> purrr::map(~(.x * t(.x)))


  list(prod = XtXs, bw = bandwidth_list)
}

#' Estimate covariance function using adaptive bandwidths
#'
#' `covariance_ll` estimates the covariance function using curves that
#' have been smoothed using an adaptive bandwidth.
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
#' - **$Gamma** Estimated covariance function
#' - **$bw** Matrix of bandwidths used to smooth each curve.
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


covariance_ll <- function(curves, grid_bandwidth =
                          lseq(0.005, 0.1, length.out = 151),
                          grid_smooth = seq(0, 1, length.out = 101),
                          k0 = 1, grid_param = seq(0.1, 0.9, length.out = 20),
                          sigma = NULL, mu0 = NULL) {

  prod_ <- smooth_curves_covariance(curves, grid_bandwidth,
                                    grid_smooth, k0,
                                    grid_param, sigma, mu0)


  wt_cond <- purrr::map(curves, ~abs(outer(.x$t, grid_smooth, "-"))) |>
    lapply(function(Ti) {
      sapply(seq_along(grid_smooth), function(t) {
        (outer(Ti[, t], prod_$bw$bw_matrix[, t], FUN = "<=") * 1) |>
          colSums(na.rm = TRUE)
      }) |>
        (\(x) (x >= k0)  * 1)()
    })

  wtws <- purrr::map(wt_cond, ~(.x * t(.x)))

  WN <- Reduce('+', wtws)

  prod_sum <- purrr::map2(wtws, prod_$prod, ~(.x * .y)) |>
    (\(x) Reduce('+', x))()

  gamma <- (1/WN) * prod_sum

  Xt_cond <- smooth_curves_mean_plugin_cov(curves, grid_smooth,
                                           k0, bandwidth = prod_$bw$bw_matrix)
  mu_ts <- purrr::map2(Xt_cond, wt_cond, ~.x * .y) |>
    (\(x) Reduce('+', x) / WN)() |>
    (\(x) x * t(x))()
  Gamma <- gamma - mu_ts

  #final step: replace diagonal band
  for (t in 1:ncol(Gamma)) {
    s <- 1
    current_cov <- Gamma[s, t - s + 1]
    while (s <= (t - s + 1)) {
      if (abs(grid_smooth[s] - grid_smooth[t - s + 1]) >
          prod_$bw$bw_matrix[s, t - s + 1]) {
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
          prod_$bw$bw_matrix[ncol(Gamma) + s - t, t]) {
        current_cov <- Gamma[ncol(Gamma) + s - t, t]
      } else {
        Gamma[ncol(Gamma) + s - t, t] <- current_cov
      }
      t <- t - 1
    }
  }
    list(Gamma = Gamma, bw = prod_$bw)
}

