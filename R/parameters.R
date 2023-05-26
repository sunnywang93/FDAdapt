#' Presmoothing with cross-validation
#'
#' Performs presmoothing of functional data, using a Nadaraya-Watson smoother
#' and a bandwidth learned from cross-validation (cv) on a subset of curves.
#' Bandwidth used is a q-th quantile of the set of those learned from cv on
#' each curve. Selected curves are uniformly sampled.
#'
#' @param data List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param N0 Number of curves randomly selected for cross-validation.
#' @param grid_smooth Grid of t2 points used to presmooth curves.
#' @param q Quantile of cv bandwidth used to presmooth curves.
#' @returns A list, containing two elements:
#' - **$presmoothed_curves** Presmoothed curves.
#' - **$cv_bandwidth** Bandwidth used.
#' @references Wang S., Patilea V., Klutchnikoff N., (2023+) - Adaptive
#' Functional Principal Components Analysis
#' @export

presmoothing_FPCA <- function(data, N0, grid_smooth, q) {

  # Randomly sample N0 curves to perform least squares cross-validation
  cv_idx <- sample(x = seq_along(data), size = N0)
  # Extract N0 curves for LSCV
  sampled_data <- data[cv_idx]
  # Perform LSCV & extract bandwidth
  lscv_h <- purrr::map_dbl(sampled_data,
                           ~np::npregbw(.x$x ~ .x$t,
                                        ckertype = "epanechnikov")$bw)
  # Select q-th quantile bandwidth
  h_cv <- unname(quantile(lscv_h, q))

  # Smooth each curve with selected bandwidth
  smoothed_curves <- smooth_curves(data = data, grid = grid_smooth, bandwidth = h_cv)

  list(presmoothed_curves = purrr::map(smoothed_curves$smoothed_curves,
                                       ~list(t = grid_smooth, x = .x)),
       cv_bandwidth = h_cv)

}


#' Smooth curves with a given bandwidth
#'
#' Performs smoothing of curves for a given bandwidth using the
#' Nadaraya-Watson smoother.
#'
#' @param data List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param grid Vector / List containing the grid of points to smooth curves. If
#' the input is a list, then smoothing is performed matching the i-th curve
#' in `data` and the i-th element on the grid.
#' @param bandwidth Numeric containing the bandwidth.
#' @returns List, containing the smoothed curves.
#' @export

smooth_curves <- function(data, grid, bandwidth) {

  if(is.list(grid)) {
    weights_num <- purrr::map2(data, grid,
                               ~epa_kernel(outer(.x$t, .y, FUN = "-") / bandwidth))
  } else {
    weights_num <- purrr::map(data,
                              ~epa_kernel(outer(.x$t, grid, FUN = "-") / bandwidth))

  }

  weights_denom <- purrr::map(weights_num, ~colSums(.x))
  weights <- purrr::map2(weights_num, weights_denom,
                            ~sweep(.x, 2, .y, FUN = "/")) |>
    rapply(f = function(x) ifelse(is.nan(x), 0, x), how = "replace")

  list(weights = weights,
       smoothed_curves = purrr::map2(weights, data, ~c(crossprod(.x, .y$x))))
}

#' Estimate moments for bandwidth constants
#'
#' Performs estimation of the second moments of the curves `Xt` and products
#' `XtXs`. These moments appear as constants in bandwidth optimization problems.
#' Moments are estimated on the same grid as the presmoothed curves, and a
#' correction is made to adjust for sparse sampling schemes.
#' @param presmoothed List, with each element containing vector of presmoothed
#' curves.
#' @returns A list, containing
#' - **$moment2** Second moment of `Xt`.
#' - **$moment2prod** Second moment of `XtXs`.
#' @references Wang S., Patilea V., Klutchnikoff N., (2023+) - Adaptive
#' Functional Principal Components Analysis
#' @export
estimate_moments <- function(presmoothed) {

  # Compute indicators for non-zero values
  indic <- purrr::map(presmoothed, ~(.x$x**2 > 0) * 1)
  sum_indic <- Reduce('+', indic)
  # Compute m2 normalised by sum of indicators
  m2 <- Reduce('+', purrr::map(presmoothed, ~.x$x**2)) / sum_indic
  m2_mean <- (Reduce('+', purrr::map(presmoothed, ~.x$x)) / sum_indic)**2
  m2_norm <- m2 - m2_mean
  # Compute indicators for c2
  indic_bi <- purrr::map(presmoothed, ~(outer(.x$x, .x$x)**2 > 0) * 1)
  sum_indic_bi <- Reduce('+', indic_bi)
  # Compute c2 normalised by sum of indicators
  c2 <- Reduce('+', purrr::map(presmoothed, ~outer(.x$x, .x$x)**2)) / sum_indic_bi
  c2_mean <- Reduce('+', purrr::map(presmoothed, ~outer(.x$x, .x$x))) / sum_indic_bi
  c2_norm <- c2 - c2_mean**2

  list(moment2 = m2_norm,
       moment2prod = c2_norm)

}

#' Estimated the local regularity based on presmoothing
#'
#' Estimates the Hurst index `H` and Hölder constant `L` using presmoothing.
#' @param data List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param grid Vector of points to estimate the regularity.
#' @param bandwidth Bandwidth used to presmooth curves.
#' @param intp Boolean, indicating whether to use interpolation method. If
#' TRUE, presmoothing is first done at the observed time points with a
#' prespecified bandwidth, before interpolating back on the t1, t2, t3 grid.
#' @references Wang S., Patilea V., Klutchnikoff N., (2023+) - Adaptive
#' Functional Principal Components Analysis
#' @export
estimate_regularity <- function(data, grid, bandwidth, intp = TRUE, gamma,
                                nknots = 15) {

  m <- purrr::map_dbl(data, ~length(.x$t)) |> mean()
  delta <- exp(-log(m)**gamma)
  t2 <- grid
  t1 <- pmax(0, t2 - delta)
  t3 <- pmin(1, t2 + delta)

  t2[t2 == 0] <- (min(t1) + min(t3)) / 2
  t2[t2 == 1] <- (max(t1) + max(t3)) / 2

  if(intp) {
    # Extract Tmi points to use as presmoothing grid plus {0, 1}
    presmooth_grid <- purrr::map(data, ~c(0, .x$t, 1))
    # Presmooth at Tmi points
    presmoothed <- smooth_curves(data = data,
                                 grid = presmooth_grid,
                                 bandwidth = bandwidth)$smoothed_curves
    # Interpolate back onto t1, t2, t3 grid
    smoothed_t1 <- purrr::map2(data, presmoothed,
                               ~interpolate1D(x = c(0, .x$t, 1),
                                              y = .y,
                                              xout = t1)$y)


    smoothed_t2 <- purrr::map2(data, presmoothed,
                               ~interpolate1D(x = c(0, .x$t, 1),
                                              y = .y,
                                              xout = t2)$y)


    smoothed_t3 <- purrr::map2(data, presmoothed,
                               ~interpolate1D(x = c(0, .x$t, 1),
                                              y = .y,
                                              xout = t3)$y)

  } else {
    # Smooth curves at t1 points
    smoothed_t1 <- smooth_curves(data, t1, bandwidth)$smoothed_curves
    # Smooth curves at t2 points
    smoothed_t2 <- smooth_curves(data, t2, bandwidth)$smoothed_curves
    # Smooth curves at t3 points
    smoothed_t3 <- smooth_curves(data, t3, bandwidth)$smoothed_curves
  }

  # Compute theta at (t1, t2)
  theta12 <- purrr::map2(smoothed_t1, smoothed_t2, ~(.x - .y)**2) |>
    (\(x) Reduce('+', x) / length(x))()

  # Compute theta at (t1, t3)
  theta13 <- purrr::map2(smoothed_t1, smoothed_t3, ~(.x - .y)**2) |>
    (\(x) Reduce('+', x) / length(x))()

  # Take the maximum between theta13 and theta12
  theta13 <- purrr::map2_dbl(theta12, theta13, ~pmax(.x, .y))

  # Compute H - Set lower and upper bounds to avoid degenerate values
  hurst <- (log(theta13) - log(theta12)) / (2 * log(2))
  hurst <- pmax(pmin(hurst, 1), 0.2)

  hurst_fit <- stats::smooth.spline(x = t2,
                                    y = hurst,
                                    nknots = nknots)

  hurst <- pmax(pmin(predict(hurst_fit, x = t2)$y, 1), 0.2)

  # Compute L2
  L_square <- theta13 / abs(t1 - t3)**(2 * hurst) |> pmax(0.1)

  list(H = hurst,
       L = sqrt(L_square))

}


#' Estimate conditional variance
#'
#' Estimates the conditional variance of curves using ordered statistics,
#' correcting for sparse regimes by selecting only points that are sufficiently
#' close enough to each other.
#'
#' @param data List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param sigma_grid Vector containing the sampling points to estimate
#' the conditional variance.
#' @param h Optional numeric containing the bandwidth h, which is only
#' used for correcting the diagonal band.
#' @param h_power Numeric, power to raise the bandwidth to when using
#' sigma for diagonal correction for the covariance function.
#' @returns Vector containing sigma.
#' @references Wang S., Patilea V., Klutchnikoff N., (2023+) - Adaptive
#' Functional Principal Components Analysis
#' @export

estimate_sigma <- function(data, sigma_grid, h, h_power = 0.9) {

  m <- purrr::map_dbl(data, ~length(.x$t)) |> mean()

  # Use a different b for the diagonal correction
  if(missing(h)) {
    b <- 0.1
  } else {
    b <- h**h_power
  }

  # Calculate the ordered indexes of time points based on absolute distance
  T_ordered <- lapply(data, function(i) {
   apply(abs(outer(i$t, sigma_grid, FUN = "-")), 2, order)
  })
  # Extract the closest and second closest indexes for each t
  T_closest_idx <- purrr::map(T_ordered, ~apply(.x, 2, function(x) x[1]))
  T_second_idx <- purrr::map(T_ordered, ~apply(.x, 2, function(x) x[2]))
  # Extract the closest and second closest time points on each curve for each t
  T_closest <- purrr::map2(data, T_closest_idx, ~.x$t[.y])
  T_second <- purrr::map2(data, T_second_idx, ~.x$t[.y])
  # Extract the closest and second closest observed points for each t
  Y_closest <- purrr::map2(data, T_closest_idx, ~.x$x[.y])
  Y_second <- purrr::map2(data, T_second_idx, ~.x$x[.y])
  #Calculate squared differences
  diff_squared <- purrr::map2(Y_closest, Y_second, ~(.x - .y)**2)

  # Compute I(t;b)
  indicators <- purrr::map(T_second, ~(abs(.x - sigma_grid) <= b) * 1)

  sum_indicators <- Reduce('+', indicators)

  # Compute sigma hat for each curve
  sigma_hat_i <- purrr::map2(diff_squared, indicators, ~.x * .y)

  # Average over all curves
  sqrt(Reduce('+', sigma_hat_i) / (2 * sum_indicators))

}


#' Estimate minimum density of sample points
#'
#' Estimate the minimum density of sampling points using a kernel-density
#' estimator, leaving out the boundary points to account for boundary bias.
#'
#' @param data List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @returns Numeric containing the minimum density.
#' @references Wang S., Patilea V., Klutchnikoff N., (2023+) - Adaptive
#' Functional Principal Components Analysis
#' @export

estimate_density <- function(data) {

  T_all <- data |> purrr::map(~.x$t) |> unlist() |> sort()
  min(density(T_all, from = 0.15, to = 0.85)$y)

}



#' Estimate parameters for FPCA
#' Performs estimate of parameters used in FPCA, in particular the local
#' regularity, moments, noise, and density of sampling points. See
#' separate internal functions for more details.
#' @param data List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param grid_points Vector containing the sampling points to estimate the
#' parameters.
#' @param n_learn Learning set used for cross-validation bandwidth in the
#' presmoothing step. See also `presmoothing_FPCA`.
#' @param h_quantile Quantile of bandwidth selected to presmooth curves.
#' @param intp Boolean, indicating whether to use interpolation method for
#' estimating Hölder parameters. If TRUE, presmoothing is first done at the
#' observed time points using a cross-validation bandwidth, before interpolating
#' back on a grid of t1, t2, t3 points used to estimate H and L.
#' @param gamma_H Numeric, power to be used in computing delta for Hurst index.
#' @param gamma_L Numeric, power to be used in computing delta for Holder constant.
#' @returns List, containing the following elements:
#' - **$t** Vector of sampling points where estimation is performed.
#' - **$H** Vector of Hurst indexes.
#' - **$L** Vector of Hölder constants.
#' - **$sigma** Vector of conditional variance / noise.
#' - **$density** Numeric containing minimum density of sampling points.
#' - **$moment2** Vector of second moments of `Xt`.
#' - **$moment2prod** Matrix of second moment of products `XtXs`.
#' @references Wang S., Patilea V., Klutchnikoff N., (2023+) - Adaptive
#' Functional Principal Components Analysis
#' @export

estimate_parameters_FPCA <- function(data, grid_points, n_learn, h_quantile,
                                     intp = TRUE, gamma_H, gamma_L) {

  # Presmooth curves
  presmoothed_curves <- presmoothing_FPCA(data = data,
                                          N0 = n_learn,
                                          grid_smooth = grid_points,
                                          q = h_quantile)
  # Estimate H and L
  hurst_estim <- estimate_regularity(data = data,
                                     grid = grid_points,
                                     bandwidth = presmoothed_curves$cv_bandwidth,
                                     intp = intp,
                                     gamma = gamma_H)$H

  const_estim <- estimate_regularity(data = data,
                                     grid = grid_points,
                                     bandwidth = presmoothed_curves$cv_bandwidth,
                                     intp = intp,
                                     gamma = gamma_L)$L

  # Estimate moments
  moments <- estimate_moments(presmoothed = presmoothed_curves$presmoothed_curves)

  # Estimate conditional variance
  noise <- estimate_sigma(data = data,
                          sigma_grid = grid_points)

  # Estimate density of time points
  sampling_density <- estimate_density(data = data)

  list(t = grid_points,
       H = hurst_estim,
       L = const_estim,
       sigma = noise,
       density = sampling_density,
       moments2 = moments$moment2,
       moments2prod = moments$moment2prod)

}



#' Estimate the variance of raw data curves, which do not require
#' smoothing.
#'
#' `estimate_variance_irreg` estimates the variance using only the
#' raw data points, possibly irregularly sampled, without pre-smoothing
#' the curves.
#'
#' @param data List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param grid_estim Vector of sampling points to estimate the variance.
#' @returns Vector, comprising the estimated variance on the grid of
#' sampling points.
#' @export
estimate_variance_irreg <- function(data, grid_estim) {

  # Find the closest point Tmi to the grid t
  idx_grid <- lapply(data, function(x) {
    sapply(grid_estim, function(t) which.min(abs(t - x$t)))
  })

  # Compute the second moment based on extracted points
  second_mom <- purrr::map2(data, idx_grid, ~.x$x[.y]**2) |>
    (\(x) Reduce('+', x) / length(x))()

  # Compute the squared mean
  mean_squared <- purrr::map2(data, idx_grid, ~.x$x[.y]) |>
    (\(x) (Reduce('+', x) / length(x))**2)()

  # Compute variance of noisy observed curves
  variance_y <- second_mom - mean_squared

  # Estimate variance of noise
  variance_noise <- estimate_sigma(data, sigma_grid = grid_estim)**2

  #Calculate the variance of true curves
  variance_y - variance_noise

}


#' Estimate the regularity of dense functional data, by only performing
#' interpolation and directly estimating H and L in the local neighbourhoods
#' defined by delta.
#'
#'
#' @param data List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param xout Vector of sampling points to estimate the regularity parameters.
#' @param method String, containing the interpolation method.
#' @param gamma_H Numeric, power to be used in computing delta for Hurst index.
#' @param gamma_L Numeric, power to be used in computinog delta for Holder
#' constant.
#' @returns List, containing the estimated H and L.
#' @export

estimate_regularity_dense <- function(data, xout, method, gamma_H, gamma_L) {

  m <- purrr::map_dbl(data, ~length(.x$t)) |> mean()

  delta_H <- exp(-log(m)**gamma_H)
  delta_L <- exp(-log(m)**gamma_L)

  t1_H <- pmax(0, xout - delta_H)
  t3_H <- pmin(1, xout + delta_H)
  t2_H <- xout
  t2_H[t2_H == 0] <- (min(t1_H) + min(t3_H)) / 2
  t2_H[t2_H == 1] <- (max(t1_H) + max(t3_H)) / 2

  t1_L <- pmax(0, xout - delta_L)
  t3_L <- pmin(1, xout + delta_L)
  t2_L <- xout
  t2_L[t2_L == 0] <- (min(t1_L) + min(t3_L)) / 2
  t2_L[t2_L == 1] <- (max(t1_L) + max(t3_L)) / 2

  intp_t2_H <- purrr::map(data, ~interpolate1D(x = .x$t,
                                               y = .x$x,
                                               xout = t2_H,
                                               type = method))


  intp_t1_H <- purrr::map(data, ~interpolate1D(x = .x$t,
                                               y = .x$x,
                                               xout = t1_H,
                                               type = method))


  intp_t3_H <- purrr::map(data, ~interpolate1D(x = .x$t,
                                               y = .x$x,
                                               xout = t3_H,
                                               type = method))

  intp_t2_L <- purrr::map(data, ~interpolate1D(x = .x$t,
                                               y = .x$x,
                                               xout = t2_L,
                                               type = method))


  intp_t1_L <- purrr::map(data, ~interpolate1D(x = .x$t,
                                               y = .x$x,
                                               xout = t1_L,
                                               type = method))


  intp_t3_L <- purrr::map(data, ~interpolate1D(x = .x$t,
                                               y = .x$x,
                                               xout = t3_L,
                                               type = method))


  theta12_H <- purrr::map2(intp_t1_H, intp_t2_H, ~(.x$y - .y$y)**2) |>
    (\(x) Reduce('+', x) / length(x))()

  theta13_H <- purrr::map2(intp_t1_H, intp_t3_H, ~(.x$y - .y$y)**2) |>
    (\(x) Reduce('+', x) / length(x))()

  # Take the maximum between theta13 and theta12
  theta13_H <- purrr::map2_dbl(theta12_H, theta13_H, ~pmax(.x, .y))

  theta13_L <- purrr::map2(intp_t1_L, intp_t3_L, ~(.x$y - .y$y)**2) |>
    (\(x) Reduce('+', x) / length(x))()

  # Compute H - Set lower and upper bounds to avoid degenerate values
  hurst <- (log(theta13_H) - log(theta12_H)) / (2 * log(2))
  hurst <- pmax(pmin(hurst, 1), 0.2)

  # Compute L2
  L_square <- theta13_L / abs(t1_L - t3_L)**(2 * hurst)

  list(H = hurst,
       L = sqrt(L_square))

}







