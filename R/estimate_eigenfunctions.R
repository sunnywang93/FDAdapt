#' Estimate adaptive bandwidths for eigenfunctions
#'
#' `estimate_bandwidth_efunctions` returns a vector of bandwidths for each
#' eigenfunction, used for smoothing curves and is optimal for
#' eigenfunction estimation (in the `L2` norm).
#'
#' @param curves List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed Points.
#' @param grid_bandwidth Grid of points to estimate bandwidth.
#' @param grid_smooth Grid of points to smooth curves.
#' @param k0 Minimum number of points curves must possess to be used in
#' estimation.
#' @param nfunctions The number of eigenfunctions to be kept.
#' @param params A tibble, containing the following parameters:
#' - **$t$** Sampling points.
#' - **$H** Estimated Hölder exponents.
#' - **$L** Estimated Hölder constants.
#' - **$sigma** Estimated noise.
#' - **$mu** Estimated density of time points.
#' estimates it.
#' @returns A vector of length nfunctions containing the adaptive bandwidth
#' for each eigenfunction.
#' @export


estimate_bandwidth_efunctions <- function(curves, grid_bandwidth, grid_smooth,
                                          k0, nfunctions, params)
  {
  cov_gkp <- covariance_norm(curves, grid_bandwidth, grid_smooth,
                             k0, params)

  eigen_elements <- normalise_eigen(cov_gkp$cov, nfunctions)

  evalues_norm <- eigen_elements$values
  efunctions_norm <- eigen_elements$vectors

  #G x H matrix
  h_alpha_t <- t(outer(grid_bandwidth, 2 * cov_gkp$constants$H, FUN = "^")) *
    (cov_gkp$constants$L**2 * cov_gkp$kernel_int)

  #G x H x J array before integrating - column recycling used
  bias_tt <- sapply(seq(nfunctions), function(j) {
    efunctions_norm[, j]**2 * h_alpha_t
  }, simplify = "array") |>
    apply(MARGIN = c(2, 3),
          FUN = function(t) pracma::trapz(grid_smooth, t))

  bias_ts <- sapply(seq(nfunctions), function(j) {
    rowSums(efunctions_norm[, -j]**2) * cov_gkp$moment2
  }) |> apply(MARGIN = 2, function(s) pracma::trapz(grid_smooth, s))

  bias_t <- sweep(bias_tt, MARGIN = 2, STATS = bias_ts, FUN = "*")

  #h_alpha_s = h_alpha_t
  bias_ss <- sapply(seq(nfunctions), function(j) {
    rowSums(efunctions_norm[, -j]**2) * h_alpha_t
  }, simplify = "array") |>
    apply(MARGIN = c(2, 3), function(s) pracma::trapz(grid_smooth, s))

  bias_st <- sapply(seq(nfunctions), function(j) {
    efunctions_norm[, j]**2 * cov_gkp$moment2
  }) |>
    apply(MARGIN = 2, function(t) pracma::trapz(grid_smooth, t))

  bias_s <- sweep(bias_ss, MARGIN = 2, STATS = bias_st, FUN = "*")

  bias_constant <- sapply(seq(nfunctions), function(j) {
    1 / sum((evalues_norm[j] - evalues_norm[-j])**2)
  }) * 2

  bias_term <- sweep(bias_s + bias_t, MARGIN = 2,
                     STATS = bias_constant, FUN = "*")

  variance_tt <- sapply(seq(nfunctions), function(j) {
    efunctions_norm[, j]**2 * cov_gkp$moment2
  })

  variance_ts <- sapply(seq(nfunctions), function(j) {
    rowSums(efunctions_norm[, -j]**2)
  })

  variance_t_num <- sapply(seq(nfunctions), function(j) {
    outer(variance_tt[, j], variance_ts[, j])
  }, simplify = "array")

  variance_t <- sapply(seq(nfunctions), function(j) {
    sapply(seq_along(grid_bandwidth), function(h) {
      variance_t_num[,,j] / cov_gkp$Ngamma_ts[,,h]
    }, simplify = "array")
  }, simplify = "array")

  variance_t[is.nan(variance_t)] <- 0

  variance_t <- variance_t |>
    apply(MARGIN = c(2, 3, 4), function(s) pracma::trapz(grid_smooth, s)) |>
    apply(MARGIN = c(2, 3), function(t) pracma::trapz(grid_smooth, t))

  #compute variance s in the same way and multiply by the common constant
  variance_s_num <- aperm(variance_t_num, perm = c(2, 1, 3))

  variance_s <- sapply(seq(nfunctions), function(j) {
    sapply(seq_along(grid_bandwidth), function(h) {
      variance_s_num[,,j] / cov_gkp$Ngamma_st[,,h]
    }, simplify = "array")
  }, simplify = "array")

  variance_s[is.nan(variance_s)] <- 0

  variance_s <- variance_s |>
    apply(MARGIN = c(2, 3, 4), function(t) pracma::trapz(grid_smooth, t)) |>
    apply(MARGIN = c(2, 3), function(s) pracma::trapz(grid_smooth, s))

  m <- purrr::map_dbl(curves, ~length(.x$t)) |> mean()

  variance_constant <- sapply(seq(nfunctions), function(j) {
    1 / sum((evalues_norm[j] - evalues_norm[-j])**2)
  }) *  max(params$sigma**2)

  variance_term <- sweep(variance_s + variance_t, MARGIN = 2,
                         STATS = variance_constant, FUN = "*")

  regularise_st <- sapply(seq_along(grid_bandwidth), function(h) {
    cov_gkp$moment2_prod * (1 / cov_gkp$WN[,,h]) - (1 / length(curves))
  }, simplify = "array")

  regularise_ts <- sapply(seq(nfunctions), function(j) {
    outer(efunctions_norm[, j]**2, rowSums(efunctions_norm[, -j]**2))
  }, simplify = "array")

  regularise <- sapply(seq(nfunctions), function(j) {
    sapply(seq_along(grid_bandwidth), function(h) {
      regularise_st[,,h] * regularise_ts[,,j]
    }, simplify = "array")
  }, simplify = "array") |>
    apply(MARGIN = c(2, 3, 4), function(s) pracma::trapz(grid_smooth, s)) |>
    apply(MARGIN = c(2, 3), function(t) pracma::trapz(grid_smooth, t))

  regularise_constant <- sapply(seq(nfunctions), function(j) {
    1 / sum((evalues_norm[j] - evalues_norm[-j])**2)
  })

  regularise_term <- sweep(regularise, MARGIN = 2,
                           STATS = regularise_constant, FUN = "*")

  risk <- bias_term + variance_term + regularise_term

  min_h_index <- apply(risk, MARGIN = 2, which.min)

  h_star <- grid_bandwidth[min_h_index]

  M_avg <- purrr::map_dbl(curves, ~length(.x$t)) |> mean()

  N_effective <- mean(cov_gkp$WN[,, min_h_index]) * M_avg

  h_constant <- log(N_effective)**(abs(log(h_star) / log(N_effective)))

  h_star * h_constant
}

#' Smooth curves with adaptive eigenfunction bandwidth
#'
#' `smooth_curves_efunctions` smooths irregularly sampled curves onto a specified
#' grid to approximate densely observed curves in functional data. Used
#' primarily to estimate eigenfunctions downstream.
#'
#' @param curves List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param grid_bandwidth Grid of points to estimate bandwidth.
#' @param grid_smooth Grid of points to smooth curves.
#' @param k0 Minimum number of points curves must possess to be used in
#' estimation.
#' @param nfunctions The number of eigenfunctions to be kept.
#' @param params A tibble, containing the following parameters:
#' - **$t$** Sampling points.
#' - **$H** Estimated Hölder exponents.
#' - **$L** Estimated Hölder constants.
#' - **$sigma** Estimated noise.
#' - **$mu** Estimated density of time points.
#' estimates it.
#' @returns A list of two elements, containing
#' * A list of curves, with
#'   - $t Sampling points.
#'   - $x Observed points.
#' * A vector of bandwidths used to smooth each curve.
#' @export


smooth_curves_efunctions <- function(curves, grid_bandwidth,
                                     grid_smooth, k0,
                                     nfunctions, params)
{
  bandwidth <- estimate_bandwidth_efunctions(curves, grid_bandwidth,
                                             grid_smooth, k0,
                                             nfunctions, params)

  M_length <- curves |> purrr::map_dbl(~(length(.x$t)))

  #dim Mi x G x J for each curve i
  Wm_num <- purrr::map(curves, ~outer(.x$t, grid_smooth, FUN = "-")) |>
    purrr::map(~outer(.x, bandwidth, FUN = "/")) |>
    purrr::map(~epa_kernel(.x))

  Wm_denom <- purrr::map(Wm_num, ~colSums(.x, na.rm = TRUE))

  Wm_num_t <- purrr::map(Wm_num, ~aperm(.x, c(2, 1, 3)))

  Wm <- purrr::map2(Wm_num_t, Wm_denom, function(x, y) {
    sapply(seq(nfunctions), function(j) {
      x[,,j] / y[, j]
    }, simplify = "array")
  })

  Wm <- lapply(Wm, function(w) replace(w, is.nan(w), 0))

  Xt <- purrr::map2(Wm, curves, function(w, y) {
    sapply(seq(nfunctions), function(j) {
      w[,,j] %*% y$x
    })
  })
  list(smoothed_curves = Xt, bandwidth = bandwidth)
}


#' Estimate eigenfunctions of functional data with adaptive bandwidths
#'
#' `efunctions_adaptive` estimates eigenfunctions using an adaptive bandwidth.
#' Curves can be irregularly sampled, possess heterogenous
#' degrees of smoothness and be relatively sparse (but should have
#' at least 20 points per curve). Uses the "smooth first, then estimate"
#' approach.
#'
#' @param curves List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param grid_bandwidth Grid of points to estimate bandwidth.
#' @param grid_smooth Grid of points to smooth curves.
#' @param k0 Minimum number of points curves must possess to be used in
#' estimation.
#' @param params A tibble, containing the following parameters:
#' - **$t$** Sampling points.
#' - **$H** Estimated Hölder exponents.
#' - **$L** Estimated Hölder constants.
#' - **$sigma** Estimated noise.
#' - **$mu** Estimated density of time points.
#' estimates it.
#' @param nfunctions The number of eigenfunctions to be kept.
#' @returns A list containing
#' * Normalised eigenfunctions.
#' * Bandwidth used for smoothing curves.
#' @references Patilea V., Wang S. (2022+) - Adaptive Functional
#' Principal Components Analysis
#' @export


efunctions_adaptive <- function(curves, grid_bandwidth, grid_smooth, k0,
                                nfunctions = 10, params)
{
  m <- purrr::map_dbl(curves, ~length(.x$t)) |> mean()

  smooth_curves <- smooth_curves_efunctions(curves, grid_bandwidth,
                                            grid_smooth, k0, nfunctions,
                                            params)

  mu_eigen <- mean_plugin_evalues(curves, smooth_curves$smoothed_curves,
                                  smooth_curves$bandwidth, grid_smooth, k0)


  centered_curves <- lapply(smooth_curves$smoothed_curves, function(i) {
    i - mu_eigen$mu
  })

  weighted_curves <- purrr::map2(centered_curves, mu_eigen$wt,
                                 ~.x * .y)


  weighted_cov <- lapply(weighted_curves, function(i) {
    sapply(seq(nfunctions), function(j) {
      i[, j] %*% t(i[, j])
    }, simplify = "array")
  }) |> (\(x) Reduce('+', x))()

  WN <- lapply(mu_eigen$wt, function(i) {
    sapply(seq(nfunctions), function(j) {
      i[, j] %*% t(i[, j])
    }, simplify = "array")
  }) |>
    (\(x) Reduce('+', x))()

  emp_cov <- weighted_cov / WN

  eelements <- lapply(seq(nfunctions), function(j) {
    normalise_eigen(emp_cov[,,j], nfunctions)
  })

  efunctions <- sapply(seq_along(eelements), function(j) {
    eelements[[j]]$vectors[, j]
  })

  list(eigenfunctions = efunctions,
       bandwidth = smooth_curves$bandwidth)
}



