#' Estimate noise level of curves
#'
#' @param data List, where each element represents a curve. Each curve
#' must be a list with two entries:
#'  * $t Sampling points.
#'  * $x Observed points.
#' @returns A number.
#' @export

estimate_sigma <- function(data) {
  Mbar <- data |> purrr::map_dbl(~length(.x$t)) |> mean()
  delta <- (1/log(Mbar)^2)
  indic <- data |>
    purrr::map(~as.numeric(abs(diff(sort(.x$t, decreasing = TRUE))) <= delta))
  diffsq <- data |> purrr::map(~diff(sort(.x$x, decreasing = TRUE))^2)
  sum_diffsq <- purrr::map2_dbl(diffsq, indic, ~sum(.x * .y))
  denom <- indic |> purrr::map_dbl(~sum(.x))
  sum_diffsq_norm <- sum_diffsq / (2 * denom)
  sqrt(mean(sum_diffsq_norm))
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
  min(density(T_all, from = 0.1, to = 0.9)$y)
}

#' Performs presmoothing of curves
#'
#' `presmoothing` performs presmoothing on irregularly sampled curves,
#' for the purpose of estimating parameters such as Hölder constants.
#' Performed using the Nadaraya-Watson estimator, with bandwidth
#' detailed in the references.
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
    sigma <- estimate_sigma(data)
  }
  if(is.null(mu0)) {
    mu0 <- estimate_density(data)
  }

  aa <- (init_b + 1) / 2 * init_b**2 * mu0
  c <- (sigma**(2*init_b) * init_L * aa**init_b)**(1 / (2*init_b + 1))
  psi_m <- (1 / m)**(init_b / (2 * init_b + 1))
  b_naive <- (c * psi_m / init_L)**(1 / init_b)

  t_list <- rbind(t1_list, t0_list, t3_list)

  theta_num <- lapply(seq(ncol(t_list)), function(t_col) {
    sapply(data, function(i) {
      (outer(i$t, t_list[, t_col], FUN = "-") / b_naive) |>
        bertin_kernel(b_naive) |>
        (\(x) (t(x) %*% i$x) / (length(i$t) * b_naive))()
      }) |> t()
    })

  theta_denom <- lapply(seq(ncol(t_list)), function(t_col) {
    sapply(data, function(i) {
      (outer(i$t, t_list[, t_col], FUN = "-") / b_naive) |>
        bertin_kernel(b_naive) |> colSums() |>
        (\(x) x / (length(i$t) * b_naive))()
      }) |> t() |> pmax(1/100)
    })

  theta <- purrr::map2(theta_num, theta_denom, ~.x / .y)

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
#' typically as output from the function `presmoothing`.
#'
#' @param presmoothed_data A list, containing
#' - **$x** Smoothed data points.
#' @param H0_list Vector containing Hölder exponents, typically as
#' output from `estimate_H0`.
#' @param M Mean number of points on each curve.
#' @returns A vector containing the estimated values at the
#' sampling points of presmoothed curves.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#' estimation of irregular mean and covariance functions.
#' @export

estimate_L0 <- function(presmoothed_data, H0_list, M) {

  V1 <- purrr::map2(presmoothed_data, H0_list,
                      ~ (.x$x[, 2] - .x$x[, 1])**2 /
                        abs(.x$t[2] - .x$t[1])**(2 * .y))

  V2 <- purrr::map2(presmoothed_data, H0_list,
                  ~ (.x$x[, 3] - .x$x[, 2])**2 /
                    abs(.x$t[3] - .x$t[2])**(2 * .y))

  V_mean <- mapply(function(x, y) (x + y) / 2, V1, V2)

  sqrt(colMeans(V_mean, na.rm = TRUE))
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
    sigma <- estimate_sigma(data)
  }
  if(is.null(mu0)) {
    mu0 <- estimate_density(data)
  }
  presmoothed <- presmoothing(data, t0_list = grid_estim, sigma = sigma,
                              mu0 = mu0)

  H0 <- estimate_H0(presmoothed)
  L0 <- estimate_L0(presmoothed, H0)

  tibble::tibble(t = grid_estim, H = H0, L = L0, sigma = sigma, mu0 = mu0)
}


