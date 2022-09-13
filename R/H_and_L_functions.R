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
                          t0_list = seq(.2, .8, l = 20),
                          init_b = 1,
                          init_L = 1,
                          sigma = NULL,
                          mu0 = NULL)
{

  NW <- function(t, T_mi, Y_mi, h, alpha = 1) {
    K <- function(x, beta){
      ((1 + beta) / (2 * beta))  * (1 - abs(x)^beta) * (abs(x) <= 1)
    }
    tt <- matrix(rep(t, length(T_mi)), ncol = length(T_mi))
    TT <- t(matrix(rep(T_mi, length(t)), ncol = length(t)))
    Kx <- 2 * K(x = (tt - TT) / h, beta = alpha)
    (Kx %*% Y_mi) / matrix(rowSums(Kx))
  }

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

  #b_naive <- (delta / round(m))**(1 / (2 * order + 1))
  #b_naive <- log(m) / m
  aa <- (init_b + 1) / 2 * init_b**2 * mu0
  c <- (sigma**(2*init_b) * init_L * aa**init_b)**(1 / (2*init_b + 1))
  psi_m <- (1 / m)**(init_b / (2 * init_b + 1))
  b_naive <- pmax(pmin((c * psi_m / init_L)**(1 / init_b), delta/4), log(m)/m)

  inner_loop <- function(i, data, t_list, b_naive, init_b) {
    sapply(data, function(x) {
      NW(t = t_list[, i],
         T_mi = x$t, Y_mi = x$x,
         h = b_naive, alpha = init_b)
    }) |> t()
  }

  t_list <- rbind(t1_list, t0_list, t3_list)
  purrr::map(1:ncol(t_list), ~list(
    t_list = t_list[,.x],
    x = inner_loop(.x, data = data, t_list = t_list,
                   b_naive = b_naive, init_b = 1)))
}

#' Performs estimation of Hölder exponent
#'
#' `estimate_H0` estimates the Hölder constant of presmoothed curves,
#' typically as output from the function `presmoothing`.
#'
#' @param data A list, containing
#' - **$x** Smoothed data points.
#' @param center TRUE/FALSE, where TRUE induces an additional
#'  centering in each term.
#' @returns A vector containing the estimated values at the
#' sampling points of presmoothed curves.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#' estimation of irregular mean and covariance functions.
#' @export

estimate_H0 <- function(data, center = FALSE){
  data |> purrr::map_dbl(function(d) {
    if(center) {
      a <- mean((d$x[, 3] - d$x[, 1] - mean(d$x[, 3] - d$x[, 1]))**2,
                na.rm = TRUE)
      b <- mean((d$x[, 2] - d$x[, 1] - mean(d$x[, 2] - d$x[, 1]))**2,
                na.rm = TRUE)
      c <- mean((d$x[, 3] - d$x[, 2] - mean(d$x[, 3] - d$x[, 2]))**2,
                na.rm = TRUE)
    }
    else {
      a <- mean((d$x[, 3] - d$x[, 1])**2, na.rm = TRUE)
      b <- mean((d$x[, 2] - d$x[, 1])**2, na.rm = TRUE)
      c <- mean((d$x[, 3] - d$x[, 2])**2, na.rm = TRUE)
    }
    max(min((2 * log(a) - log(b * c)) / log(16), 1), 0.1)
  }
  )
}


#' Performs estimation of Hölder constant
#'
#' `estimate_L0` estimates the Hölder constant of presmoothed curves,
#' typically as output from the function `presmoothing`.
#'
#' @param data A list, containing
#' - **$x** Smoothed data points.
#' @param H0_list Vector containing Hölder exponents, typically as
#' output from `estimate_H0`.
#' @param M Mean number of points on each curve. Defaults to NULL,
#' which calculates it.
#' @param center TRUE/FALSE, where TRUE induces an additional
#'  centering in each term.
#' @returns A vector containing the estimated values at the
#' sampling points of presmoothed curves.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#' estimation of irregular mean and covariance functions.
#' @export

estimate_L0 <- function(data, H0_list, M = NULL, center = FALSE) {
  if(is.null(M)) {
    M <- data |> purrr::map_dbl(~length(.x$x)) |> mean()
  }
  H0 <- H0_list |> purrr::map_dbl(~ .x - 1 / log(M)**1.01)
  if(center) {
    V1 <- data |>
      purrr::map2(H0,
                  ~ (.x$x[, 2] - .x$x[, 1] - mean(.x$x[, 2] - .x$x[, 1]))**2 /
                    abs(.x$t[2] - .x$t[1] - mean(.x$t[2] - .x$t[1]))**(2 * .y))
    V2 <- data |>
      purrr::map2(H0,
                  ~ (.x$x[, 3] - .x$x[, 2] - mean(.x$x[, 3] - .x$x[, 2]))**2 /
                    abs(.x$t[3] - .x$t[2] - mean(.x$t[3] - .x$t[2]))**(2 * .y))
  }
  else {
    V1 <- data |>
      purrr::map2(H0,
                  ~ (.x$x[, 2] - .x$x[, 1])**2 / abs(.x$t[2] - .x$t[1])**(2 * .y))
    V2 <- data |>
      purrr::map2(H0,
                  ~ (.x$x[, 3] - .x$x[, 2])**2 / abs(.x$t[3] - .x$t[2])**(2 * .y))
  }
  V_mean <- V1 |> purrr::map2_dfc(V2, ~ (.x + .y) / 2)
  unname(sqrt(colMeans(V_mean, na.rm = TRUE)))
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
  presmoothed <- data |>
    presmoothing(t0_list = grid_estim, sigma = sigma, mu0 = mu0)

  H0 <- estimate_H0(presmoothed)
  L0 <- estimate_L0(presmoothed, H0)

  presmoothed_2nd <- data |>
    presmoothing(t0_list = grid_estim, init_b = H0, init_L = L0,
                 sigma = sigma, mu0 = mu0)

  H1 <- estimate_H0(presmoothed_2nd)
  L1 <- estimate_L0(presmoothed_2nd, H1)
  tibble::tibble(t = grid_estim, H = H1, L = L1, sigma = sigma, mu0 = mu0)
}


