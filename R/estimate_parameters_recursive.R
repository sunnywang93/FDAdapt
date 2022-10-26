#' Performs estimation of parameters
#'
#' `estimate_holder_quantities` estimates the parameters used for
#' downstream analysis, typically for functions such as
#' `evalues_adaptive`. Cross validation to determine the local
#'  bandwidth for presmoothing, before a correction is performed
#'  based on theoretical developments.
#'
#' @param curves A list containing the raw data points with two entries:
#' - **$t** Sampling points.
#' - **$x** Observed Points
#' @param grid_param Grid of points to estimate parameters.
#' @returns A list with two elements:
#' - **$params** Tibble of parameters.
#' - **$smoothed_curves** Smoothed_curves.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#' estimation of irregular mean and covariance functions.
#' @export

estimate_holder_quantities <- function(curves, grid_param) {

  sigma <- estimate_sigma(curves, grid_param)
  #replace with np density estimator
  mu0 <- estimate_density(curves)

  m <- mean(purrr::map_dbl(curves, ~length(.x$t)))

  delta <- min(log(m)**(-1.1) / 2, 0.1)

  t2_list <- grid_param
  t1_list <- grid_param - delta
  t3_list <- grid_param + delta

  t_list <- rbind(t1_list, t2_list, t3_list)

  smoothed_curves_list <- locpoly_smooth(curves, t_list)
  smoothed_curves <- smoothed_curves_list$presmoothed
  H0 <- estimate_H(smoothed_curves)
  L0 <- estimate_L(smoothed_curves, t_list, H0)

  tibble::tibble(t = grid_param, H = H0, L = L0,
                 sigma = sigma, mu0 = mu0)
}


#' Performs presmoothing for the purpose of parameter estimation
#'
#' `locpoly_smooth` epresmooths raw curves, for the purpose of
#' parameter estimation. Cross-validation is first done to determine
#' a local bandwidth, before a correction is performed based on thereotical
#' grounding.
#'
#' @param curves A list containing the raw data points with two entries:
#' - **$t** Sampling points.
#' - **$x** Observed Points
#' @param grid Grid of points to estimate parameters. Must be a matrix, with
#' the rows representing the 3 time points, `t1, t2, t3`, and the columns
#' are the points on the estimation grid.
#' @returns An `G x N x 3` array, where `G` is the number of points on the
#' the grid and `N` is the number of curves.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#' estimation of irregular mean and covariance functions.
#' @export

#grid: t1, t2, t3 as a matrix
locpoly_smooth <- function(curves, grid) {

  sorted_grid <- sort(as.vector(grid), index.return = TRUE)

  bandwidth_list <- purrr::map(curves,
                               ~lokern::lokerns(x = .x$t, y = .x$x,
                                        x.out = sorted_grid$x)$bandwidth)

  bandwidth_list_bounded <- purrr::map(bandwidth_list,
                                       ~ifelse(.x <= 0.01, 0.01, .x))

  bandwidth_t2 <- purrr::map(bandwidth_list_bounded,
                             ~.x[sorted_grid$ix %in%
                                   seq(2, length(sorted_grid$x), by = 3)])

  bandwidth_t2_rep <- purrr::map(bandwidth_t2, ~rep(.x, each = 3))


  presmoothed <- lapply(seq_along(curves), function(i) {
    lokern::lokerns(curves[[i]]$t, curves[[i]]$x, x.out = sorted_grid$x,
                    inputb = TRUE,
                    bandwidth = bandwidth_t2_rep[[i]])
  })


  presmoothed_values <- purrr::map(presmoothed,
                                   ~list(t = .x$x.out,
                                         x = .x$est))

  presmoothed_t1 <- sapply(presmoothed_values,
                           function(x) x$x[sorted_grid$ix %in%
                                             seq(1, length(sorted_grid$x),
                                                 by = 3)])

  presmoothed_t2 <- sapply(presmoothed_values,
                           function(x) x$x[sorted_grid$ix %in%
                                             seq(2, length(sorted_grid$x),
                                                 by = 3)])

  presmoothed_t3 <- sapply(presmoothed_values,
                           function(x) x$x[sorted_grid$ix %in%
                                             seq(3, length(sorted_grid$x),
                                                 by = 3)])

  list(presmoothed = abind::abind(presmoothed_t1, presmoothed_t2,
                                  presmoothed_t3, along = 3),
       bandwidth = bandwidth_t2_rep)
}


epa_smoother <- function(curves, grid_points, bandwidth) {

  weights <- lapply(curves, function(i) {
    epa_kernel(outer(i$t, grid_points, FUN = "-") / bandwidth)
  })

  weights_denom <- purrr::map(weights, ~colSums(.x))

  #normalise and define 0 / 0 = 0
  weights_norm <- purrr::map2(weights, weights_denom,
                              ~sweep(.x, MARGIN = 2, STATS = .y, FUN = "/")) |>
    rapply(f = function(x) ifelse(is.nan(x), 0, x), how = "replace")

  mapply(function(x, y) t(x) %*% y$x, weights_norm, curves)
}

estimate_theta <- function(curves_u, curves_v) {
  mean((curves_u - curves_v)**2)
}


#' Estimated the local regularity based on presmoothing
#'
#' Estimates the local regularity with supplied presmoothed curves.
#'
#' @param smoothed_curves An array with dimensions `G x N x 3`, where
#' the rows are the number of points on the grid used for estimation,
#' the columns are the curves, and the 3rd dimension represent the 3 time points.
#' Should use `locpoly_smooth` as an auxiliary function to obtain the proper
#' inputs.
#' @returns A vector, whose length is equal to the number of sampling
#' points the presmoothed curves takes values in.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#' estimation of irregular mean and covariance functions.
#' @export
estimate_H <- function(smoothed_curves) {
  # theta1 <- apply(smoothed_curves, MARGIN = 3,
  #                 function(x) estimate_theta(x[1, ], x[3, ]))
  #
  # theta2 <- apply(smoothed_curves, MARGIN = 3,
  #                 function(x) estimate_theta(x[1, ], x[2, ]))
  #
  # H <- (log(theta1) - log(theta2)) / 2*log(2)
  # pmin(pmax(H, 0.1), 1)
  theta1 <- rowMeans((smoothed_curves[,,1] - smoothed_curves[,,3])**2)
  theta2 <- rowMeans((smoothed_curves[,,2] - smoothed_curves[,,3])**2)
  H <- (log(theta1) - log(theta2)) / 2 * log(2)
  pmin(pmax(H, 0.1), 1)
}

#' Estimated the local Hölder constant based on presmoothing
#'
#' Estimates the local Hölder constant with supplied presmoothed curves.
#'
#' @param smoothed_curves An array with dimensions `G x N x 3`, where
#' the rows are the number of points on the grid used for estimation,
#' the columns are the curves, and the 3rd dimension represent the 3 time points.
#' Should use `locpoly_smooth` as an auxiliary function to obtain the proper
#' inputs.
#' @param t_list Must be a matrix, with
#' the rows representing the 3 time points, `t1, t2, t3`, and the columns
#' are the points on the estimation grid.
#' @param H A vector with estimated regularity from `estimate_H`, with
#' the same length as the number of columns in `t_list`.
#' @returns A vector, whose length is equal to the number of sampling
#' points the presmoothed curves takes values in.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#' estimation of irregular mean and covariance functions.
#' @export
estimate_L <- function(smoothed_curves, t_list, H) {

  theta1 <- rowMeans((smoothed_curves[,,2] - smoothed_curves[,,3])**2)
  theta2 <- rowMeans((smoothed_curves[,,2] - smoothed_curves[,,1])**2)

  L1 <- theta1 / abs(t_list[3, ] - t_list[2, ])**(2 * H)
  L2 <- theta2 / abs(t_list[2, ] - t_list[1, ])**(2 * H)

  sqrt((L1 + L2) / 2)
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
#' - **sigma** Estimated noise.
#' - **mu0** Estimated density of sampling points.
#' @param grid_smooth Grid of points to smooth curves.
#' @returns A list, containing
#' - **$VarXt** Estimated Variance of `Xt`.
#' - **$VarXtXs** Estimated Variance of `XtXs`.
#' - **$EXt2** Estimated second moment of `Xt`.
#' - **$EXtXs2** Estimated second momenet of `XtXs`.
#' @export
estimate_variance_curves <- function(data, params, grid_smooth) {

  sigma <- params$sigma

  mu0 <- params$mu0

  m <- data |> purrr::map_dbl(~length(.x$t)) |> mean()

  grid_tibble <- params

  X_hat <- sapply(data, function(i) lokern::lokerns(x = i$t,
                                                    y = i$x,
                                                    x.out = grid_smooth)$est)

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


