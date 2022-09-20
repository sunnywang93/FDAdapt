#' Estimate noise level of curves using only observed points
#'
#' #' `estimate_sigma` estimates the variance of curves, using only
#' information from curves. Used as a preliminary estimate of sigma,
#' usually as a first pass to the function `estimate_sigma_recursive`.
#'
#' @param data List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @returns A number.
#' @export

estimate_sigma <- function(data) {
  Mbar <- data |> purrr::map_dbl(~length(.x$t)) |> mean()
  delta <- log(log(Mbar)) / Mbar
  indic <- data |>
    purrr::map(~as.numeric(abs(diff(sort(.x$t, decreasing = TRUE))) <= delta))
  diffsq <- data |> purrr::map(~diff(sort(.x$x, decreasing = TRUE))^2)
  sum_diffsq <- purrr::map2_dbl(diffsq, indic, ~sum(.x * .y))
  denom <- indic |> purrr::map_dbl(~sum(.x))
  sum_diffsq_norm <- sum_diffsq / (2 * denom)
  sqrt(mean(sum_diffsq_norm))
}

#' Estimate noise level of curves with presmoothing
#'
#' `estimate_sigma_recursive` estimates the variance of curves, using
#' a presmoothing approach. Presmoothing is done using Bertin's approach,
#' detailed in the references. Has a minimum bandwidth corresponding to
#' 5% of the interval if the calculated bandwidth is too big.
#'
#' @param data List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @returns A number.
#' #' @references Bertin L, (2004) - Minimax exact constant in sup-norm for
#' nonparametric regression with random design.
#' @export

estimate_sigma_recursive <- function(data) {
  Mbar <- data |> purrr::map_dbl(~length(.x$t)) |> mean()
  interval <- purrr::map_dbl(data, ~(max(.x$t) - min(.x$t))) |> max()
  bandwidth_min <- interval * 0.075
  sigma <- estimate_sigma(data)
  mu0 <- estimate_density(data)
  bandwidth <- bertin_bandwidth(sigma, mu0, init_b = 1, init_L = 1, m = Mbar) |>
    pmin(bandwidth_min)
  presmooth_curves <- bertin_smoother_recursive(data, bandwidth)
  resid <- purrr::map2(data, presmooth_curves, ~(.x$x - .y$x)**2)
  purrr::map(resid, ~mean(.x))
  unlist(resid) |> mean(na.rm = TRUE) |> sqrt()
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
                          mu0 = NULL)
{
  m <- data |> purrr::map_dbl(~ length(.x$t)) |> mean()

  delta <- min(log(m)^(-1.1), 0.2)
  t1_list <- t0_list - delta / 2
  t3_list <- t0_list + delta / 2

  if(is.null(sigma)) {
    sigma <- estimate_sigma_recursive(data)
  }
  if(is.null(mu0)) {
    mu0 <- estimate_density(data)
  }

  # interval_length <- purrr::map_dbl(data, ~max(.x$t) - min(.x$t)) |> max()
  # bandwidth_min <- interval_length * 0.05

  bandwidth <- bertin_bandwidth(sigma, mu0, init_b, init_L, m) #|>
    #pmin(bandwidth_min) |> pmax(log(m) / m)

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

  H0_smooth <- H0_list |> purrr::map_dbl(~.x - 1/log(M)**1.01)
  #H0_smooth <- stats::smooth.spline(x = t0_list, y = H0_list, nknots = 8)$y

  V1 <- purrr::map2_dbl(presmoothed_data, H0_smooth,
                      ~ mean((.x$x[, 2] - .x$x[, 1])**2, na.rm = TRUE) /
                        abs(.x$t[2] - .x$t[1])**(2 * .y))

  V2 <- purrr::map2_dbl(presmoothed_data, H0_smooth,
                  ~ mean((.x$x[, 3] - .x$x[, 2])**2, na.rm = TRUE) /
                    abs(.x$t[3] - .x$t[2])**(2 * .y))

  sqrt((V1 + V2) / 2)
}


#' Performs twice recursive estimation of parameters
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
#' @returns Tibble with estimated parameters as columns.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#' estimation of irregular mean and covariance functions.
#' @export

estimate_holder_const <- function(data,
                                  grid_estim = seq(0.2, 0.8, length.out = 20),
                                  sigma = NULL,
                                  mu0 = NULL) {

  if(is.null(sigma)) {
    sigma <- estimate_sigma_recursive(data)
  }
  if(is.null(mu0)) {
    mu0 <- estimate_density(data)
  }

  presmoothed <- presmoothing(data, t0_list = grid_estim, sigma = sigma,
                              mu0 = mu0)

  H0 <- estimate_H0(presmoothed)

  M <- purrr::map_dbl(data, ~length(.x$t)) |> mean()

  L0 <- estimate_L0(presmoothed, H0, grid_estim, M)

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
#' - **$EXtXs2** Estimated second momenet of `XtXs`.
#' @export
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
#' using `estimate_sigma` or `estimate_sigma_recursive`.
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


#' Smooths curves on its own grid
#'
#' `bertin_smoother_recursive` performs curve smoothing using a kernel,
#' referred to as the Bertin kernel (see `bertin_kernel`). Used only
#' as an auxiliary function to `estimate_sigma_recursive`. Smoothing
#' is done on the grid of each curve's observed points.
#'
#' @param data List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @param bandwidth Bandwidth used to smooth curves, usually
#'  calculated from `bertin_bandwidth`.
#' @returns A list, containing
#' - **$t** Sample points of curves.
#' - **$x** Smoothed data points.
#' @references Bertin L, (2004) - Minimax exact constant in sup-norm for
#' nonparametric regression with random design.
#' @export

bertin_smoother_recursive <- function(data, bandwidth) {

  theta_num <- lapply(data, function(i) {
    (outer(i$t, i$t, FUN = "-") / bandwidth) |>
      bertin_kernel(bandwidth) |>
      (\(x) (x %*% i$x) / (length(i$t) * bandwidth))() |>
      c()
  })

  theta_denom <- lapply(data, function(i) {
    (outer(i$t, i$t, FUN = "-") / bandwidth) |>
      bertin_kernel(bandwidth) |> colSums() |>
      (\(x) x / (length(i$t) * bandwidth))() |>
      pmax(1/100)
  })

  theta <- purrr::map2(theta_num, theta_denom, ~(.x / .y))

  purrr::map2(data, theta, ~list(t = .x$t,
                                 x = .y))

}




