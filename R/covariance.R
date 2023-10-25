#' Computes the optimal bandwidth for the covariance in terms of MISE
#'
#' Given a list of functional data, the optimal bandwidth for the covariance
#' function in terms of mean integrated squared error (MISE) is computed on
#' a grid of bandwidths. A plug-in bandwidth rule based on risk bounds is
#' computed.
#'
#' @param curves List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param params List, containing
#' - **$t** Vector of evaluation points on which the parameters were estimated on.
#' - **$H** Vector containing the estimated Hölder exponent.
#' - **$L** Vector containing the estimated Hölder constant.
#' - **$sigma** Vector containing the estimated noise.
#' - **$moments2** Vector containing the estimated variance of the curves.
#' - **$moments2prod** Matrix containing the variance of the product of curves.
#' @param h_grid Vector containing the grid of bandwidths.
#' @param t_grid Vector containing the evaluation points on which the curves will
#' be smoothed on.
#' @param inflate Boolean, indicating whether to inflate the bandwidths by
#' the discretization errors.
#' @returns Numeric, containing the optimal bandwidth.
#' @export

bw_cov <- function(curves, params, h_grid, t_grid, inflate = TRUE) {

  # Compute constant kernel
  cst_kernel <- 1.5 * (1 / (1 + 2*params$H) - 1 / (3 + 2*params$H))

  # Compute bias term
  bias_gamma <- 4 * bias_kernel(H = params$H,
                                L = params$L,
                                cst = cst_kernel,
                                h_grid = h_grid) |>
    apply(1, function(x) pracma::trapz(t_grid, x)) *
    pracma::trapz(t_grid, params$moments2)

  # Compute variance term
  Ngamma_list <- N_gamma(curves = curves, x = t_grid, bw = h_grid)

  var_num <- outer(params$moments2, params$sigma**2)

  var_gamma <- sapply(seq_along(h_grid), function(h) {
    var_num / Ngamma_list$Ngamma[,,h]
  }, simplify = "array")

  var_gamma <- apply(var_gamma,
                     c(2, 3),
                     function(x) pracma::trapz(t_grid, x)) |>
    apply(2, function(x) pracma::trapz(t_grid, x))

  var_gamma2 <- sapply(seq_along(h_grid), function(h) {
    var_num / t(Ngamma_list$Ngamma[,,h])
  }, simplify = "array")

  var_gamma2 <- apply(var_gamma2,
                      c(2, 3),
                      function(x) pracma::trapz(t_grid, x)) |>
    apply(2, function(x) pracma::trapz(t_grid, x))

  variance_gamma <- 2 * var_gamma + 2 * var_gamma2

  # Compute regularising term
  reg_gamma <- sapply(seq_along(h_grid), function(h) {
    params$moments2prod *
      ((1 / Ngamma_list$WN_bi[,,h]) - (1 / length(curves)))
  }, simplify = "array") |>
    apply(c(2, 3), function(x) pracma::trapz(t_grid, x)) |>
    apply(2, function(x) pracma::trapz(t_grid, x))

  # Compute final risk
  risk_gamma <- bias_gamma + variance_gamma + reg_gamma

  idx_min <- which.min(risk_gamma)

  # Find and return optimal bandwidth
  if(inflate) {
    log(1 / h_grid[idx_min]) * h_grid[idx_min]
  } else {
    h_grid[idx_min]
  }

}

#' Computes the covariance function based on MISE
#'
#' Given a list of functional data, the covariance function is computed based on
#' a plug-in bandwidth rule which minimises the mean-integrated squared error (MISE).
#' An adaptive approach is taken by first estimating the regularity of the sample
#' paths, and then computing the risk bounds based on the estimated regularity.
#'
#' @param Y_list List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param xout Vector containing the evaluation points for smoothing curves.
#' @param hout Vector containing the grid of bandwidths to optimise over.
#' @param param_list List of estimated parameters, for example from the
#' function `estimate_parameters_FPCA`.
#' @param h_power Numeric, power to raise bandwidth to when estimating sigma
#' for the purposes of diagonal correction for the covariance function.
#' @param inflate_bw Boolean, indicating whether to inflate the bandwidths by
#' the discretization errors.
#' @param diag_cor Boolean, indicating whether to perform diagonal correction.
#' @returns Matrix containing the estimated covariance function.
#' @export

cov_ISE <- function(Y_list, xout, hout, param_list,
                    h_power = 0.9, inflate_bw = TRUE, diag_cor = TRUE) {

  if(length(param_list$t) != length(xout)) {
    stop("param_list must be on the same grid as grid_smooth!")
  }


  bw_gamma <- bw_cov(curves = Y_list,
                     params = param_list,
                     h_grid = hout,
                     t_grid = xout,
                     inflate = inflate_bw)


  wi <- curves_select(curves = Y_list,
                      x = xout,
                      h = bw_gamma)

  WN_bi <- apply(wi, 2, function(x) tcrossprod(x), simplify = FALSE) |>
    (\(x) Reduce('+', x))()

  X_hat <- smooth_curves(data = Y_list,
                         grid = xout,
                         bandwidth = bw_gamma)

  mu_hat <- purrr::imap(X_hat$smoothed_curves, ~wi[, .y] * .x) |>
    (\(x) Reduce('+', x))() / rowSums(wi)

  Gamma_hat <- purrr::imap(X_hat$smoothed_curves,
                           ~tcrossprod(.x - mu_hat) * tcrossprod(wi[, .y])) |>
    (\(x) Reduce('+', x))() / WN_bi

  if(diag_cor) {

    d_hat <- diagonal_bias(curves = Y_list,
                           t_grid = xout,
                           h = bw_gamma,
                           zeta = h_power,
                           wm_list = X_hat$weights,
                           W = wi,
                           WN = WN_bi)

    Gamma_hat <- Gamma_hat - d_hat

  }


  list(
    bw = bw_gamma,
    gamma = Gamma_hat
  )

}










