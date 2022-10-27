#' Estimate adaptive bandwidths for eigenvalues
#'
#' `estimate_bandwidth_evalues` returns a vector of bandwidths for each
#' eigenvalue, used for smoothing curves and is optimal for eigenvalue estimation.
#'
#' @param curves List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed Points.
#' @param grid_bandwidth Grid of points to estimate bandwidth.
#' @param grid_smooth Grid of points to smooth curves.
#' @param k0 Minimum number of points curves must possess to be used in
#' estimation.
#' @param grid_param Grid of points to estimate parameters.
#' @param sigma Noise level of curves.
#' @param mu0 Density lower bound for time points.
#' @param nvalues The number of eigenvalues to be kept.
#' @returns A vector of length nvalues containing the adaptive bandwidth
#' for each eigenvalue.
#' @examples
#' estimate_bandwidth_evalues(curves = curves_list,
#' grid_bandwidth = seq(0, 1, length.out = 151),
#' grid_smooth = seq(0, 1, length.out = 101),
#' k0 = 1,
#' grid_param = seq(0, 1, length.out = 20),
#' sigma = 0.1, mu0 = 1, nvalues = 10)
#' @export

estimate_bandwidth_evalues <- function(curves, grid_bandwidth, grid_smooth, k0,
                                       nvalues, params)
  {

  cov_norm <- covariance_norm(curves, grid_bandwidth, grid_smooth,
                             k0, params)

  #obtain the normalised eigenvalues
  eigen_elements <- normalise_eigen(cov_norm$cov, nvalues)

  #extract normalised eigenvalues
  evalues_norm <- eigen_elements$values
  efunctions_norm <- eigen_elements$vectors

  h_alpha <- sapply(grid_bandwidth, function(h) {
    cov_norm$constants$L**2 * h**(2*cov_norm$constants$H)
  })

  #G x H matrix for fixed j - column recycling used here
  bias_t <- sapply(seq(nvalues), function(j) {
    efunctions_norm[, j]**2 * cov_norm$kernel_int * h_alpha
  }, simplify = "array")

  #integrate over the grid: output H x J
  bias_t_int <- apply(bias_t,
                      MARGIN = c(2, 3),
                      FUN = function(t) pracma::trapz(grid_smooth, t))


  bias_s <- sapply(seq(nvalues), function(j) {
    cov_norm$moment2 * efunctions_norm[, j]**2
  })

  #output: J x 1
  bias_s_int <- apply(bias_s, MARGIN = 2,
                      function(t) pracma::trapz(grid_smooth, t))

  #output: J x H
  #equal to s and t, with 2 factor in front of each
  #so we must multiply by 4 at the end
  bias_term <- apply(bias_t_int, MARGIN = 1,
                     function(h) h * bias_s_int)

  #loop over h and j - column recycling over m2
  variance_s <- sapply(seq_along(grid_bandwidth), function(h) {
    sapply(seq(nvalues), function(j) {
      outer(efunctions_norm[, j]**2, efunctions_norm[, j]**2) *
        cov_norm$moment2 / cov_norm$Ngamma_st[,,h]
    }, simplify = "array")
  }, simplify = "array")

  variance_s[is.nan(variance_s)] <- 0

  variance_t <- aperm(variance_s, c(2, 1, 3, 4))
  #perform integration iteratively
  variance_s_int <- apply(variance_s,
                          MARGIN = c(2, 3, 4),
                          function(var) pracma::trapz(grid_smooth, var)) |>
    apply(MARGIN = c(2, 3), function(var) pracma::trapz(grid_smooth, var))

  rm(variance_s)

  variance_t_int <- apply(variance_t,
                          MARGIN = c(2, 3, 4),
                          function(var) pracma::trapz(grid_smooth, var)) |>
    apply(MARGIN = c(2, 3), function(var) pracma::trapz(grid_smooth, var))

  rm(variance_t)

  m <- purrr::map_dbl(curves, ~length(.x$t)) |> mean()

  variance_term <- (max(sigma**2) * variance_s_int +
    max(sigma**2) * variance_t_int)

  regularising_term <- sapply(seq_along(grid_bandwidth), function(h) {
    sapply(seq(nvalues), function(j) {
      outer(efunctions_norm[, j]**2, efunctions_norm[, j]**2) *
        cov_norm$moment2_prod * (1 / cov_norm$WN[,, h] - 1 / length(curves))
    }, simplify = "array")
  }, simplify = "array") |>
    apply(MARGIN = c(2, 3, 4), function(t) pracma::trapz(grid_smooth, t)) |>
    apply(MARGIN = c(2, 3), function(s) pracma::trapz(grid_smooth, s))

  risk <- 4 * bias_term + variance_term + regularising_term

  min_h_index <- apply(risk, MARGIN = 1, which.min)

  h_star <- sapply(min_h_index, function(id) grid_bandwidth[id])

  M_avg <- purrr::map_dbl(curves, ~length(.x$t)) |> mean()

  N_effective <- mean(cov_norm$WN[,, min_h_index]) * M_avg

  h_constant <- log(N_effective)**(abs(log(h_star) / log(N_effective)))

  h_star * h_constant
}

#' Smooth curves with adaptive eigenvalue bandwidth
#'
#' `smooth_curves_evalues` smooths irregularly sampled curves onto a specified
#' grid to approximate densely observed curves in functional data. Used
#' primarily to estimate eigenvalues downstream.
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
#' @param sigma Noise level of curves.
#' @param mu0 Density lower bound for time points.
#' @param nvalues The number of eigenvalues to be kept.
#' @returns A list of two elements, containing
#' * A list of curves, with
#'   - $t Sampling points.
#'   - $x Observed points.
#' * A vector of bandwidths used to smooth each curve.
#' @examples
#' smooth_curves_evalues(curves = curves_list,
#' grid_bandwidth = seq(0, 1, length.out = 151),
#' grid_smooth = seq(0, 1, length.out = 101),
#' k0 = 1,
#' grid_param = seq(0, 1, length.out = 20),
#' sigma = 0.1, mu0 = 1, nvalues = 10)
#' @export

smooth_curves_evalues <- function(curves, grid_bandwidth, grid_smooth, k0,
                                  nvalues, params)
{

  bw_vector <- estimate_bandwidth_evalues(curves, grid_bandwidth, grid_smooth, k0,
                                          nvalues, params)

  M_length <- curves |> purrr::map_dbl(~(length(.x$t)))

  #dim Mi x G x J for each curve i
  Wm_num <- purrr::map(curves, ~outer(.x$t, grid_smooth, FUN = "-")) |>
    purrr::map(~outer(.x, bw_vector, FUN = "/")) |>
    purrr::map(~epa_kernel(.x))

  Wm_denom <- purrr::map(Wm_num, ~colSums(.x, na.rm = TRUE))

  Wm_num_t <- purrr::map(Wm_num, ~aperm(.x, c(2, 1, 3)))

  Wm <- purrr::map2(Wm_num_t, Wm_denom, function(x, y) {
    sapply(seq(nvalues), function(j) {
      x[,,j] / y[, j]
    }, simplify = "array")
  })

  Wm <- lapply(Wm, function(w) replace(w, is.nan(w), 0))

  Xt <- purrr::map2(Wm, curves, function(w, y) {
    sapply(seq(nvalues), function(j) {
      w[,,j] %*% y$x
    })
  })
  list(smoothed_curves = Xt, bw = bw_vector)
}


#' Estimate eigenvalues of functional data with adaptive bandwidths
#'
#' `evalues_adaptive` estimates eigenvalues using an adaptive bandwidth.
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
#' @param grid_param Grid of points to estimate parameters.
#' @param sigma Noise level, if known. Defaults to NULL, which estimates it.
#' @param mu0 Density lower bound for time points. Defaults to NULL, which
#' estimates it.
#' @param nvalues The number of eigenvalues to be kept.
#' @returns A list containing
#' * Normalised eigenvalues.
#' * Bandwidth used for smoothing curves.
#' @examples
#' evalues_adaptive(curves = curves_list,
#' grid_bandwidth = seq(0, 1, length.out = 151),
#' grid_smooth = seq(0, 1, length.out = 101),
#' k0 = 1,
#' grid_param = seq(0, 1, length.out = 20),
#' sigma = 0.1, mu0 = 1, nvalues = 10)
#' @references Patilea V., Wang S. (2022+) - Adaptive Functional
#' Principal Components Analysis
#' @export

evalues_adaptive <- function(curves, grid_bandwidth, grid_smooth, k0,
                             nvalues = 10, params) {

  smooth_curves <- smooth_curves_evalues(curves, grid_bandwidth,
                                         grid_smooth, k0, nvalues,
                                         params)

  mu_eigen <- mean_plugin_evalues(curves, smooth_curves$smoothed_curves,
                                  smooth_curves$bw, grid_smooth, k0)


  centered_curves <- lapply(smooth_curves$smoothed_curves, function(i) {
    i - mu_eigen$mu
  })

  weighted_curves <- purrr::map2(centered_curves, mu_eigen$wt,
                                 ~.x * .y)


  weighted_cov <- lapply(weighted_curves, function(i) {
    sapply(seq(nvalues), function(j) {
      i[, j] %*% t(i[, j])
    }, simplify = "array")
  }) |> (\(x) Reduce('+', x))()

  WN <- lapply(mu_eigen$wt, function(i) {
    sapply(seq(nvalues), function(j) {
      i[, j] %*% t(i[, j])
    }, simplify = "array")
  }) |>
    (\(x) Reduce('+', x))()

  emp_cov <- weighted_cov / WN

  eelements <- lapply(seq(nvalues), function(j) {
    normalise_eigen(emp_cov[,,j], nvalues)
  })

  evalues <- sapply(seq_along(eelements), function(j) {
    eelements[[j]]$values[j]
  })

  list(eigenvalues = evalues,
       bandwidth = smooth_curves$bw)
}



