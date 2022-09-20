#' Generate exponential grid for optimisation
#'
#' @param from Starting point of grid.
#' @param to End point of grid.
#' @param length.out Number of points in the grid.
#' @returns A vector containing grid values.
#' @examples
#' lseq(from = 1, to = 100, length.out = 51)
#' @export

lseq <- function (from = 1, to = 100, length.out = 51) {
  exp(seq(log(from), log(to), length.out = length.out))
}

#' @export
epa_kernel <- function(y) {
  3/4 * (1 - y^2) * (abs(y) <= 1)
}


#' @export
bertin_kernel <- function(x, beta){
  ((1 + beta)/ (2 * beta))  * (1 - abs(x)^beta) * (abs(x) <= 1)
}


#' Estimate adaptive bandwidth for mean estimation
#'
#' `estimate_bandwidth_mean` estimates the adaptive bandwidth
#' used for mean function estimation.
#'
#' @param data List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param params Tibble of parameters, containing
#' - **$H** Estimated regularity.
#' - **$L** Estimated Hölder constant.
#' - **$sigma** Noise level of curves.
#' - **$mu0** Density lower bound for time points. Defaults to NULL,
#' which estimates it.
#' @param grid_bandwidth Grid of points to estimate bandwidth.
#' @param grid_smooth Grid of points to smooth curves.
#' @param k0 Minimum number of points curves must possess to be used
#' in estimation.
#' @returns A tibble, containing the input `params` and vector of
#' bandwidths.
#' @export

estimate_bandwidth_mean <- function(curves, params, grid_bandwidth,
                                    grid_smooth,
                                    k0) {


  H_grid_smooth <- pracma::interp1(x = c(0, params$t, 1),
                                   y = c(params$H[1], params$H,
                                         params$H[length(params$H)]),
                                   xi = grid_smooth, method = "linear")

  L_grid_smooth <- pracma::interp1(x = c(0, params$t, 1),
                                   y = c(params$L[1], params$L,
                                         params$L[length(params$L)]),
                                   xi = grid_smooth, method = "linear")

  grid_tibble <- tibble::tibble(t = grid_smooth,
                                H = H_grid_smooth,
                                L = L_grid_smooth,
                                sigma = max(params$sigma))

  #calculate the constants needed to smooth bandwidth
  #here we use the one associated with epanechnikov kernel
  cst_kernel <- 1.5 * (1/(1 + 2 * grid_tibble$H) - 1/(3 + 2 * grid_tibble$H))
  #cst_kernel <- 1/4 * gamma(grid_tibble$H + 1) + 1/4 * gamma(grid_tibble$H + 2)
  q1 <- grid_tibble$L/factorial(floor(grid_tibble$H)) * sqrt(cst_kernel)
  q2 <- grid_tibble$sigma
  moments <- estimate_variance_curves(data = curves,
                                      params = params,
                                      grid_smooth = grid_smooth,
                                      sigma = grid_tibble$sigma,
                                      mu0 = max(params$mu0))
  q3 <- sqrt(moments$varXt)

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
  wi <- sapply(Mi_vec, function(Mi) {
    indic[Mi,,] |> colSums(na.rm = TRUE) |>
      (\(x) x >= k0)() * 1
  }, simplify = "array") |>
    aperm(c(1, 3, 2))

  #dim G x H
  WN <- wi |> aperm(c(2, 1, 3)) |> colSums(na.rm = TRUE)

  #compute max Wm
  #only compute over max of m
  #dim G x N X H
  max_Wm_num <- sapply(Mi_vec, function(Mi) {
    T_diff[Mi, ] |> abs() |> matrixStats::colMins()
  }, simplify = "array") |> outer(grid_bandwidth, FUN = "/") |>
    epa_kernel()

  K_sum <- T_diff |>
    outer(grid_bandwidth, FUN = "/") |>
    epa_kernel()

  Wm_denom <- sapply(Mi_vec, function(Mi) {
    K_sum[Mi, ,] |> colSums(na.rm = TRUE)
  }, simplify = "array") |> aperm(c(1, 3, 2))

  Wm <- max_Wm_num / Wm_denom
  Wm[is.nan(Wm)] <- 0

  #each matrix has dim (G x N)
  Ni <- wi / Wm
  Ni[is.nan(Ni)] <- 0

  #dim G x H
  WN_inv2 <- 1/WN^2

  wiNi <- wi / Ni
  wiNi[is.nan(wiNi)] <- 0
  #dim G x H
  wiNi_sum <- wiNi |> aperm(c(2, 1, 3)) |> colSums(na.rm = TRUE)

  Nmu <- 1/(WN_inv2*wiNi_sum)
  #Nmu[is.nan(Nmu)] <- 0

  #column recycling - dim G x H for all q terms
  q2_term <- q2^2 / Nmu
  #q2_term[!is.finite(q2_term)] <- 0

  q1_term <- sapply(grid_bandwidth, function(h) {
    q1^2 * h^(2*grid_tibble$H)
  })
  #column recycling
  q3_term <- q3^2  * (1/WN - 1/length(curves))

  risk <- q1_term + q2_term + q3_term

  min_h_index <- apply(risk, 1, which.min)

  grid_tibble |> dplyr::mutate(bandwidth = grid_bandwidth[min_h_index])
}



#' Smooth curves using adaptive bandwidth for mean function estimation
#'
#' `smooth_curves_mean` smooths each curve using an input bandwidth,
#' typically the one returned by `estimate_bandwidth_mean`. Smoothing is
#' targeted for the purpose of mean function estimation.
#'
#' @param curves List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param params Tibble of parameters, containing
#' - **$H** Estimated regularity.
#' - **$L** Estimated Hölder constant.
#' - **$sigma** Noise level of curves.
#' - **$mu0** Density lower bound for time points. Defaults to NULL,
#' which estimates it.
#' @param grid_bandwidth Grid of points to estimate bandwidth.
#' @param grid_smooth Grid of points to smooth curves.
#' @param k0 Minimum number of points curves must possess to be used
#' in estimation.
#' @returns A list, containing
#' - **$params** Tibble with estimated parameters.
#' - **$x_hat** Matrix of curves, with columns indexing each curve and
#' rows indexing each point on `grid_smooth`.
#' @export

smooth_curves_mean <- function(curves, params, grid_bandwidth,
                               grid_smooth, k0) {

  params <- estimate_bandwidth_mean(curves = curves, params = params,
                                    grid_bandwidth = grid_bandwidth,
                                    grid_smooth = grid_smooth,
                                    k0 = k0)

  M_length <- curves |> purrr::map_dbl(~(length(.x$t)))
  M_length_cum <- cumsum(M_length)

  vec <- seq(1, max(M_length_cum))
  vec_cut <- cut(vec, breaks = c(1, M_length_cum), labels = FALSE,
                 include.lowest = TRUE)
  Mi_vec <- split(vec, vec_cut) |> unname()

  T_diff <- lapply(curves, function(x) x$t) |>
    unlist() |>
    outer(grid_smooth, FUN = "-")

  indic <- T_diff |>
    abs() |>
    (\(x) x <= matrix(params$bandwidth, nrow = sum(M_length),
                      ncol = length(grid_smooth), byrow = TRUE))() * 1

  #dim G x N
  wi <- sapply(Mi_vec, function(Mi) {
    indic[Mi,] |> colSums(na.rm = TRUE) |>
      (\(x) x >= k0)() * 1
  })

  #dim G x 1
  WN <- wi |> rowSums(na.rm = TRUE)

  Wm_num <- (T_diff / matrix(params$bandwidth, nrow = sum(M_length),
                        ncol = length(grid_smooth), byrow = TRUE)) |>
    epa_kernel()

  #transposing makes the next calculation easier due to alignment of
  #column indices - N x G
  Wm_denom <- sapply(Mi_vec, function(Mi) {
    Wm_num[Mi, ] |> colSums(na.rm = TRUE)
  }) |> t()

  #we repeat the matrix to have the same dimensions as the numerator
  #for easy element by element division
  Wm_denom_rep <- matrix(data = rep(NA, times = sum(M_length) * length(grid_smooth)),
                         nrow = sum(M_length),
                         ncol = length(grid_smooth))

  for(i in seq_along(M_length)) {
    Wm_denom_rep[Mi_vec[[i]], ] <- matrix(Wm_denom[i, ],
                                   nrow = M_length[i],
                                   ncol = length(grid_smooth),
                                   byrow = TRUE)
  }
  #dim M x G
  Wm <- Wm_num / Wm_denom_rep
  Wm[is.nan(Wm)] <- 0

  Wm_reshaped <- lapply(Mi_vec, function(Mi) {
    Wm[Mi, ]
  })
  #dim G x N
  Xt <- sapply(seq_along(curves), function(i) {
    t(Wm_reshaped[[i]]) %*% curves[[i]]$x
  })

  #final output: N x G matrix
  list(params = params, x_hat = t(Xt))
}

#' Estimate the mean function using adaptive bandwidths
#'
#' `mean_ll` estimates the mean from curves smoothed using an
#' adaptive bandwidth, typically the one returned by `smooth_curves_mean`.
#'
#' @param data List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param grid_bandwidth Grid of points to estimate bandwidth.
#' @param grid_smooth Grid of points to smooth curves.
#' @param k0 Minimum number of points curves must possess to be used
#' in estimation.
#' @returns A list, containing
#' - **$params** Tibble with estimated parameters.
#' - **$mu_hat** Vector containing the mean function at each point
#' `t` in `grid_smooth`.
#' @export

mean_ll <- function(data, grid_bandwidth = lseq(0.005,
                                                0.1, length.out = 151),
                    grid_smooth = seq(0, 1, length.out = 101),
                    k0 = 2) {

  params <- estimate_holder_const(curves, sigma = sigma, mu0 = mu0,
                                  grid_estim = grid_param)

  smoothed <- smooth_curves_mean(curves = data, params = params,
                                 grid_bandwidth = grid_bandwidth,
                                 grid_smooth = grid_smooth, k0 = k0)
  mu_hat <- smoothed$x_hat |> colMeans(na.rm = TRUE)
  list(params = smoothed$params, mu_hat = mu_hat)
}

#' Smooth curves using a plugin bandwidth for covariance estimation
#'
#' `smooth_curves_mean_plugin_cov` smooths each curve using a plugin
#' bandwidth, typically to estimate the mean used to center the
#' estimated covariance function. Used as an auxiliary function
#' in `covariance_ll`.
#'
#' @param curves List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param grid_smooth Grid of points to smooth curves.
#' @param k0 Minimum number of points curves must possess to be used
#' in estimation.
#' @param bandwidth Matrix of bandwidths for each point `(s,t)`.
#' @returns A list containing smoothed curves,
#' @export

smooth_curves_mean_plugin_cov <- function(curves, grid_smooth,
                                      k0, bandwidth) {

  M_length <- curves |> purrr::map_dbl(~(length(.x$t)))

  #dim Mi x G
  Wm_num <- purrr::map(curves, ~outer(.x$t, grid_smooth, "-")) |>
    lapply(function(Ti) {
      sapply(seq_along(grid_smooth), function(t) {
        outer(Ti[, t], bandwidth[, t], FUN = "/")
      }, simplify = "array") |> epa_kernel()
    })

  Wm_denom <- purrr::map(Wm_num, ~colSums(.x, na.rm = TRUE)) |>
    purrr::map2(M_length, ~replicate(.y, .x)) |>
    purrr::map(~aperm(.x, c(3, 1, 2)))

  Wm <- purrr::map2(Wm_num, Wm_denom, ~(.x / .y) |>
                      (\(x) replace(x, is.nan(x), 0))())

  Xt <- purrr::map2(curves, Wm, function(i, w) {
    sapply(seq_along(grid_smooth), function(s) {
      t(w[,,s]) %*% i$x
    })
  })

  return(Xt)
}

#' Computes the mean function using a plugin bandwidth for
#' eigenvalue estimation
#'
#' `mean_plugin_evalues` computes the mean function used to center
#' the smoothed curves for the purpose of eigenvalue estimation.
#'
#' @param curves List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param smoothed_curves List containing the smoothed curves, with
#' each element containing the values of one curve.
#' @param bandwidth Matrix of bandwidths for each point `(s,t)`.
#' @param grid_smooth Grid of points to smooth curves.
#' @param k0 Minimum number of points curves must possess to be used
#' in estimation.
#' @returns A list containing
#' - **$mu** Mean curve.
#' - **$wt** Weights `wi` used to compute the mean function.
#' @export

mean_plugin_evalues <- function(curves, smoothed_curves, bandwidth,
                                grid_smooth, k0) {

  indic <- lapply(curves, function(i) {
    outer(i$t, grid_smooth, FUN = "-") |> abs() |>
      outer(bandwidth, FUN = "<=") * 1
  })

  wt <- lapply(indic, function(i) {
    colSums(i, na.rm = TRUE) |>
      (\(x) (x >= k0) * 1)()
  })

  WN <- Reduce('+', wt)

  sum_curves <- purrr::map2(wt, smoothed_curves, ~(.x * .y)) |>
    (\(x) Reduce('+', x))()

  list(mu = sum_curves / WN,
       wt = wt)
}






