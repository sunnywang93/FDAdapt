#' Estimate pointwise noise level of curves using only observed points
#'
#' `estimate_sigma` estimates the variance of curves, using only
#' information from curves.
#'
#' @param data List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param grid_param Vector containing the sampling points to estimate
#' the noise.
#' @returns Vector with the same length as `grid_param`.
#' @export

estimate_sigma <- function(data, grid_param = seq(.1, .9, length.out = 20)) {
  Mbar <- data |> purrr::map_dbl(~length(.x$t)) |> mean()
  delta <- 2 * sqrt(log(Mbar)) / Mbar
  #diffsq <- data |> purrr::map(~diff(sort(.x$x, decreasing = TRUE))**2)
  idx <- data |> purrr::map(~order(.x$t))
  diffsq <- data |> purrr::map2(idx, ~c(0, diff(.x$x[.y])**2))
  #indic <- purrr::map2(data, idx, ~(abs(diff(.x$t[.y])) <= delta) * 1)
  indic <- purrr::map(data,
                      ~(abs(outer(.x$t, grid_param, FUN = "-")) <= delta) * 1)
  sum_diffsq <- purrr::map2(diffsq, indic, ~t(.y) %*% .x)
  denom <- purrr::map(indic, ~colSums(.x, na.rm = TRUE))
  sum_diffsq_norm <- purrr::map2(sum_diffsq, denom, ~.x / (2 * .y))
  modifiedSum <- function(x, y) {
    replace(x, is.na(x), 0) + replace(y, is.na(y), 0)
  }
  c(sqrt(Reduce(modifiedSum, sum_diffsq_norm) / length(sum_diffsq_norm)))
}


#' Estimate minimum density of sample points
#'
#' `estimate_density` estimates the minimum density of time points using
#' the kernel density estimator.
#'
#' @param data List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @returns A number.
#' @export

estimate_density <- function(data) {
  T_all <- data |> purrr::map(~.x$t) |> unlist() |> sort()
  min(density(T_all, from = 0.15, to = 0.85)$y)
}

#' Performs presmoothing of curves
#'
#' `presmoothing` performs presmoothing on irregularly sampled curves,
#' for the purpose of estimating parameters such as Hölder constants.
#' Performed using a modified Nadaraya-Watson estimator, with bandwidth
#' detailed in the references. A lower and upper bound on the bandwidths
#' are imposed in order to avoid degenerate cases.
#'
#' @param data List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param t0_list Vector of sampling points which presmoothing
#'  is performed.
#' @param init_b Initialised Hölder exponent.
#' @param init_L Initialised Hölder constant.
#' @param sigma Noise level, if known. Defaults to NULL, which estimates it.
#' @param mu0 Density lower bound for time points. Defaults to NULL,
#'  which estimates it.
#' @returns A list, containing
#' - **$t_list** Sample points used for presmoothing.
#' - **$x** Smoothed data points.
#' @references Bertin L, (2004) - Minimax exact constant in sup-norm for
#' nonparametric regression with random design.
#' @export

presmoothing <- function (data,
                          t0_list = seq(.1, .9, l = 20),
                          init_b = 1,
                          init_L = 1,
                          sigma = NULL,
                          mu0 = NULL,
                          k0_list = 2)
{
  m <- data |> purrr::map_dbl(~ length(.x$t)) |> mean()

  delta <- min(log(m)^(-2.5), 0.2)
  t1_list <- t0_list - delta / 2
  t3_list <- t0_list + delta / 2

  if(is.null(sigma)) {
    sigma <- estimate_sigma(data)
  }
  if(is.null(mu0)) {
    mu0 <- estimate_density(data)
  }

  interval_length <- purrr::map_dbl(data, ~max(.x$t) - min(.x$t)) |> max()
  bandwidth_min <- interval_length * 0.15

  # H0_order <- estimate_H0_order_list(data, t0_list = t0_list,
  #                                    k0_list = k0_list, sigma = sigma)

  # bandwidth <- bertin_bandwidth(sigma, mu0, init_b = H0_order, init_L, m) |>
  #   mean()

  bandwidth <- bertin_bandwidth(sigma, mu0, init_b, init_L, m)|>
    pmin(bandwidth_min)

  t_list <- rbind(t1_list, t0_list, t3_list)

  theta <- lapply(seq(ncol(t_list)), function(t_col) {
    bertin_smoother(data, t_list[, t_col], bandwidth)
  })

  purrr::map(1:ncol(t_list), ~list(t_list = t_list[, .x],
                                   x = theta[[.x]]))
}

#' Performs estimation of Hölder exponent
#'
#' `estimate_H0` estimates the Hölder constant of presmoothed curves,
#' typically as output from the function `presmoothing`.
#'
#' @param presmoothed_data A list, containing
#' - **$x** Smoothed data points.
#' @returns A vector containing the estimated values at the
#' sampling points of presmoothed curves.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#' estimation of irregular mean and covariance functions.
#' @export

estimate_H0 <- function(presmoothed_data){
  presmoothed_data |> purrr::map_dbl(function(d) {
    a <- mean((d$x[, 3] - d$x[, 1])**2, na.rm = TRUE)
    b <- mean((d$x[, 2] - d$x[, 1])**2, na.rm = TRUE)
    c <- mean((d$x[, 3] - d$x[, 2])**2, na.rm = TRUE)

    max(min((2 * log(a) - log(b * c)) / log(16), 1), 0.1)
  }
  )
}


#' Performs estimation of Hölder constant
#'
#' `estimate_L0` estimates the Hölder constant of presmoothed curves,
#' typically as output from the function `presmoothing`. Smooths the
#' input values of `H0` by splines to get more stable results.
#'
#' @param presmoothed_data A list, containing
#' - **$x** Smoothed data points.
#' @param H0_list Vector containing Hölder exponents, typically as
#' output from `estimate_H0`.
#' @param t0_list Vector containing the sampling points of `H0_list`.
#' @param M Mean number of points on each curve.
#' @returns A vector containing the estimated values at the
#' sampling points of presmoothed curves.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#' estimation of irregular mean and covariance functions.
#' @export

estimate_L0 <- function(presmoothed_data, H0_list, t0_list, M) {

  #H0_smooth <- H0_list |> purrr::map_dbl(~.x - 1/log(M)**1.01)
  #H0_smooth <- stats::smooth.spline(x = t0_list, y = H0_list, nknots = 8)$y
  H0_smooth <- H0_list

  V1 <- purrr::map2_dbl(presmoothed_data, H0_smooth,
                      ~ mean((.x$x[, 2] - .x$x[, 1])**2, na.rm = TRUE) /
                        abs(.x$t[2] - .x$t[1])**(2 * .y))

  V2 <- purrr::map2_dbl(presmoothed_data, H0_smooth,
                  ~ mean((.x$x[, 3] - .x$x[, 2])**2, na.rm = TRUE) /
                    abs(.x$t[3] - .x$t[2])**(2 * .y))

  sqrt((V1 + V2) / 2)
}


#' Estimates local regularity using splines
#'
#' Estimates the local regularity `H(t)` by presmoothing, where
#' presmoothing is performed using splines on irregularly sampled curves.
#'
#' @param data List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param t0_list Vector of sampling points which presmoothing
#'  is performed.
#' @returns A vector containing `H(t)` evaluated at `t0_list`.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#' estimation of irregular mean and covariance functions.
#' @export


presmoothing_splines <- function(data, t0_list = seq(.1, .9, l = 20)) {

  m <- data |> purrr::map_dbl(~length(.x$t)) |> mean()
  delta <- min(log(m)^(-1.1), 0.2)

  t1_list <- t0_list - delta / 2
  t3_list <- t0_list + delta / 2
  t_list <- c(t1_list, t0_list, t3_list)

  spline_model <- purrr::map(data,
                             ~stats::smooth.spline(x = .x$t, y = .x$x,
                                                   cv = TRUE))

  y_pred <- purrr::map(spline_model, ~predict(.x, x = sort(t_list)))

  idx1 <- seq(1, length(t_list), by = 3)
  idx2 <- seq(2, length(t_list), by = 3)
  idx3 <- seq(3, length(t_list), by = 3)

  V1 <- purrr::map(y_pred, ~(.x$y[idx1] - .x$y[idx3])**2) |>
    (\(x) Reduce('+', x) / length(x))()

  V2 <- purrr::map(y_pred, ~(.x$y[idx1] - .x$y[idx2])**2) |>
    (\(x) Reduce('+', x) / length(x))()

  H0 <- pmax(pmin((log(V1) - log(V2)) / 2*log(2), 1), 0.1)

  Q1 <- V1 / abs(t1_list - t3_list)**(2 * H0)

  Q2 <- V2 / abs(t0_list - t3_list)**(2 * H0)

  L0 <- (sqrt(Q1) + sqrt(Q2)) / 2

  tibble::tibble(H = H0, L = L0)
}

#' Estimates regularity at a point
#'
#' Estimates the regularity of a curve at a fixed point `t0`.
#' Usually used as an auxiliary function for `estimate_H0_order_list` to
#' estimate the local regularity on a grid of points, and should be used
#' only for curves that have regularity less than 1.
#'
#' @param data A list containing the raw data points with two entries:
#' - **$t** Sampling points.
#' - **$x** Observed Points
#' @param t0 Time point at which `H(t0)` should be computed.
#' @param k0 Number of closest neighbours to be used in estimation.
#' If `k0 < 2`, the regularity will be set to 1.
#' @param sigma Noise level of the curves, assumed to be known. If
#' unknown, one can use the `estimate_sigma` function to estimate it.
#' @return Numeric for `H(t0)`.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2022) - Learning the
#' smoothness of noisy curves with application to online curve estimation.
#' @export


estimate_H0_order <- function(data, t0 = 0, k0 = 2, sigma) {

  num_first <- 2 * log(2)
  num_second <- 0
  denom <- 2 * log(2)

  #obtain index of 4k0 - 3th closest point to t0 on each curve
  idx_4k <- purrr::map_dbl(data, ~order(abs(.x$t - t0))[4 * k0 - 3])

  #obtain index of 2k0 - 1th closest point to t0 on each curve
  idx_2k <- purrr::map_dbl(data, ~order(abs(.x$t - t0))[2 * k0 - 1])

  #obtain index of k0th closest point to t0 on each curve
  idx_k <- purrr::map_dbl(data, ~order(abs(.x$t - t0))[k0])

  #compute 2k-1th order statistic
  a <- purrr::imap_dbl(data, ~(.x$x[idx_4k[.y]] - .x$x[idx_2k[.y]])**2) |>
    mean(na.rm = TRUE)

  #compute 4k-3th order statistic
  b <- purrr::imap_dbl(data, ~(.x$x[idx_2k[.y]] - .x$x[idx_k[.y]])**2) |>
    mean(na.rm = TRUE)

  if ((a > 2 * sigma**2) & (b > 2 * sigma**2) & (a > b)) {
    num_first <- log(a - 2 * sigma**2)
    num_second <- log(b - 2 * sigma**2)
  }

  (num_first - num_second) / denom
}


#' Estimates regularity on a grid of points
#'
#' Estimates the local regularity on a grid of points
#' `t0_list`, using order statistics (see references). Makes use
#' of the function `estimate_H0_order` which estimates
#' the local regularity at a point `t0`. Should only be used for curves
#' with a regularity less than one.
#'
#' @param data A list containing the raw data points with two entries:
#' - **$t** Sampling points.
#' - **$x** Observed Points
#' @param t0_list A vector containing the sampling points at which
#' `H(t)` should be estimated.
#' @param k0_list The number of neighbours to be considered in
#' estimation corresponding to the vector `t0_list`. To set
#' the same number of neighbours used in estimation, use a numeric.
#' @param sigma Noise level of the curves, assumed to be known. If
#' unknown, one can use the `estimate_sigma` function to estimate it.
#' @return A vector of points.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2022) - Learning the
#' smoothness of noisy curves with application to online curve estimation.


estimate_H0_order_list <- function(data, t0_list, k0_list, sigma) {
  purrr::map2_dbl(t0_list, k0_list,
                  ~estimate_H0_order(data, t0 = .x, k0 = .y, sigma = sigma))
}


#' Estimates Hölder constant at a fixed point
#'
#' Estimates the Hölder constant at a fixed point `t0`
#' using order statistics (see references). Should only be used for curves
#' with a regularity less than one.
#'
#' @param data A list containing the raw data points with two entries:
#' - **$t** Sampling points.
#' - **$x** Observed Points
#' @param t0 Numeric, the sampling point to estimate at which `L(t0)`
#' should be estimated.
#' @param H0 Numeric, the regularity `H(t0)`.
#' @param k0 Number of neighbours to be used in estimation.
#' @param sigma Noise level of the curves, assumed to be known. If
#' unknown, one can use the `estimate_sigma` function to estimate it.
#' @return A numeric.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2022) - Learning the
#' smoothness of noisy curves with application to online curve estimation.

estimate_L0_order <- function(data, t0, H0, k0, sigma) {
  num <- 1
  denom <- 1

  #obtain index of 4k0 - 3th closest point to t0 on each curve
  idx_4k <- purrr::map_dbl(data, ~order(abs(.x$t - t0))[4 * k0 - 3])

  #obtain index of 2k0 - 1th closest point to t0 on each curve
  idx_2k <- purrr::map_dbl(data, ~order(abs(.x$t - t0))[2 * k0 - 1])

  #obtain index of k0th closest point to t0 on each curve
  idx_k <- purrr::map_dbl(data, ~order(abs(.x$t - t0))[k0])

  #compute 2k-1th order statistic
  a <- purrr::imap_dbl(data, ~(.x$x[idx_4k[.y]] - .x$x[idx_2k[.y]])**2) |>
    mean(na.rm = TRUE)

  #compute 4k-3th order statistic
  b <- purrr::imap_dbl(data, ~(.x$x[idx_2k[.y]] - .x$x[idx_k[.y]])**2) |>
    mean(na.rm = TRUE)

  c <- purrr::imap_dbl(data,
                       ~abs(.x$t[idx_4k[.y]] - .x$t[idx_2k[.y]])**(2 * H0)) |>
    mean(na.rm = TRUE)

  d <- purrr::imap_dbl(data,
                       ~abs(.x$t[idx_2k[.y]] - .x$t[idx_k[.y]])**(2 * H0)) |>
    mean(na.rm = TRUE)

  if((a > b) & (c > d)) {
    num <- a - b
    denom <- c - d
  }
  sqrt(num / denom)
}

#' Estimates Hölder constant on a grid of points
#'
#' Estimates the Hölder constant on a grid of points `t0_list`
#' using order statistics (see references). Makes use of the function
#' `estimate_L0_order`. Should only be used for curves
#' with a regularity less than one.
#'
#' @param data A list containing the raw data points with two entries:
#' - **$t** Sampling points.
#' - **$x** Observed Points
#' @param t0_list Vector, containing the sampling points at which
#' `L(t)` should be estimated.
#' @param H0_list Vector, containing the points `H(t)` estimated on the
#' grid `t0_list`. One can use `estimate_H0_order_list` to estimate `H(t)`
#' at `t0_list`.
#' @param k0_list The number of neighbours to be considered in estimation
#' corresponding to the vector t0_list. To set the same number of
#' neighbours used in estimation, use a numeric.
#' @param sigma Noise level of the curves, assumed to be known. If
#' unknown, one can use the `estimate_sigma` function to estimate it.
#' @return A vector of points.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2022) - Learning the
#' smoothness of noisy curves with application to online curve estimation.


estimate_L0_order_list <- function(data, t0_list, k0_list, H0_list, sigma) {
  purrr::pmap_dbl(list(t0_list, H0_list, k0_list), function(t0, H0, k0)
    estimate_L0_order(data, t0 = t0, k0 = k0, H0 = H0, sigma = sigma))
}



#' Performs estimation of parameters
#'
#' `estimate_holder_const` estimates the parameters used for
#' downstream analysis, typically for functions such as
#' `evalues_adaptive`.
#'
#' @param data A list containing the raw data points with two entries:
#' - **$t** Sampling points.
#' - **$x** Observed Points
#' @param grid_estim Grid of points to estimate parameters.
#' @param sigma Noise level, if known. Defaults to NULL, which
#' estimates it.
#' @param mu0 Density lower bound for time points. Defaults to NULL,
#' which estimates it.
#' @param presmoothing Presmoothing option. Can be either "splines",
#' "bertin" or FALSE. If FALSE, uses the order statistics approach instead
#' (see `estimate_H0_order`).
#' @returns Tibble with estimated parameters as columns.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#' estimation of irregular mean and covariance functions.
#' @export

estimate_holder_const <- function(data,
                                  grid_estim = seq(0.2, 0.8, length.out = 20),
                                  sigma = NULL,
                                  mu0 = NULL,
                                  presmoothing = "splines") {

  if(is.null(sigma)) {
    sigma <- estimate_sigma(data)
  }
  if(is.null(mu0)) {
    mu0 <- estimate_density(data)
  }

  if(presmoothing == "splines") {
    constants <- presmoothing_splines(data, t0_list = grid_estim)
    H0 <- constants$H
    L0 <- constants$L
  }

  else if(presmoothing == "bertin") {
    presmoothed <- presmoothing(data, t0_list = grid_estim, sigma = sigma,
                                mu0 = mu0)

    H0 <- estimate_H0(presmoothed)

    M <- purrr::map_dbl(data, ~length(.x$t)) |> mean()

    L0 <- estimate_L0(presmoothed, H0, grid_estim, M)
  }

  else {
    H0 <- estimate_H0_order_list(data, t0_list = grid_estim,
                                 k0_list = 2, sigma = sigma)

    L0 <- estimate_L0_order_list(data, t0_list = grid_estim,
                                 k0_list = 2, H0_list = H0,
                                 sigma = sigma)
  }
  tibble::tibble(t = grid_estim, H = H0, L = L0, sigma = sigma, mu0 = mu0)
}


#' Estimate moments of curves for downstream analysis
#'
#' `estimate_variance_curves` estimates the moments used in
#' downstream adaptive estimation.
#'
#' @param data List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param params Tibble of parameters, containing
#' - **$H** Estimated regularity.
#' - **$L** Estimated Hölder constant.
#' @param grid_smooth Grid of points to smooth curves.
#' @returns A list, containing
#' - **$VarXt** Estimated Variance of `Xt`.
#' - **$VarXtXs** Estimated Variance of `XtXs`.
#' - **$EXt2** Estimated second moment of `Xt`.
estimate_variance_curves <- function(data, params, grid_smooth,
                                     sigma = NULL, mu0 = NULL) {
  if(is.null(sigma)) {
    sigma <- estimate_sigma(data)
  }

  if(is.null(mu0)) {
    mu0 <- estimate_density(data)
  }

  m <- data |> purrr::map_dbl(~length(.x$t)) |> mean()

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
                                sigma = max(sigma))


  #now smooth each curve using Bertin's bandwidth
  bandwidth <- bertin_bandwidth(sigma, mu0, init_b = 1, init_L = 1, m)
  X_hat <- bertin_smoother(data, grid_smooth, bandwidth) |> t()

  E_Xt2 <- apply(X_hat**2, 1, mean, na.rm = TRUE)
  var_Xt <- E_Xt2 - (apply(X_hat, 1, mean, na.rm = TRUE)**2)

  #G X G matrix for each curve
  X_hat_prod <- lapply(seq_along(data), function(i) {
    X_hat[, i] %*% t(X_hat[, i])
  })

  #use modifiedSum to ensure NAs do not affect the computation
  modifiedSum <- function(x, y) {
    replace(x, is.na(x), 0) + replace(y, is.na(y), 0)
  }

  X_bar_prod <- Reduce(modifiedSum, X_hat_prod) / length(X_hat_prod)

  diff_var_prod <- lapply(X_hat_prod, function(i) {
    (i - X_bar_prod)**2
  })

  var_XtXs <- Reduce(modifiedSum, diff_var_prod)/length(diff_var_prod)

  EXtXs2 <- lapply(X_hat_prod, function(i) {
    i^2
  }) |> (\(x) Reduce(modifiedSum, x) / length(x))()

  list(varXt = var_Xt, varXtXs = var_XtXs, EXt2 = E_Xt2,
       EXtXs2 = EXtXs2)
}

#' Calculates bandwidth for smoothing
#'
#' `bertin_bandwidth` calculates the bandwidth used for downstream smoothing.
#' Approach is detailed in the references.
#'
#' @param sigma Noise level of the curves. If unknown, can be estimated from
#' using `estimate_sigma`.
#' @param mu0 Minimum density of time points. If unknown, can be estimated
#' using `estimate_density`.
#' @param init_b Initialised Hölder exponent.
#' @param init_L Initialised Hölder constant.
#' @param m Average number of sampling points per curve.
#' @returns A scalar or vector, depending on the inputs.
#' @references Bertin L, (2004) - Minimax exact constant in sup-norm for
#' nonparametric regression with random design.
#' @export

bertin_bandwidth <- function(sigma, mu0, init_b, init_L, m) {
  aa <- (init_b + 1) / 2 * init_b**2 * mu0
  c <- (sigma**(2*init_b) * init_L * aa**init_b)**(1 / (2*init_b + 1))
  psi_m <- (1 / m)**(init_b / (2 * init_b + 1))
  (c * psi_m / init_L)**(1 / init_b)
}

#' Smooths curves using modified Nadaraya-Watson smoother
#'
#' `bertin_smoother` performs curve smoothing using a kernel, referred to as
#' the Bertin kernel (see `bertin_kernel`). Usually used as an auxiliary
#' function for others such as `presmoothing`.
#'
#' @param data List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param grid Vector of sampling points which smoothing
#'  is performed.
#' @param bandwidth Bandwidth used to smooth curves, usually
#'  calculated from `bertin_bandwidth`.
#' @returns A matrix, with the rows representing curves
#' and columns the time points where smoothing is performed.
#' @references Bertin L, (2004) - Minimax exact constant in sup-norm for
#' nonparametric regression with random design.
#' @export

bertin_smoother <- function(data, grid, bandwidth) {

  theta_num <- sapply(data, function(i) {
      (outer(i$t, grid, FUN = "-") / bandwidth) |>
        bertin_kernel(bandwidth) |>
        (\(x) (t(x) %*% i$x) / (length(i$t) * bandwidth))()
    }) |> t()

  theta_denom <- sapply(data, function(i) {
      (outer(i$t, grid, FUN = "-") / bandwidth) |>
        bertin_kernel(bandwidth) |> colSums() |>
        (\(x) x / (length(i$t) * bandwidth))()
    }) |> t()

  theta <- theta_num / theta_denom
  theta[is.nan(theta)] <- 0
  theta
}




