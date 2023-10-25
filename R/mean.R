
#' Computes the optimal bandwidth for the mean function
#'
#' The optimal bandwidth in terms of mean integrated squared error (MISE)
#' for the purpose of mean function estimation is computed.
#'
#' @param curves List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param params List of parameters, which can be obtained from
#' `estimate_parameters_FPCA`.
#' @param h_grid Vector of points for bandwidth optimization.
#' @param t_grid Vector of sampling points to smooth curves.
#' @returns Numeric, containing the optimal bandwidth for the mean function.
#' @export

bw_mean <- function(curves, param_smooth, h_grid, t_grid) {


  # Compute constant kernel
  cst_kernel <- 1.5 * (1 / (1 + 2*param_smooth$H) - 1 / (3 + 2*param_smooth$H))

  # Compute bias term
  bias_mu <- bias_kernel(H = param_smooth$H,
                         L = param_smooth$L,
                         cst = cst_kernel,
                         h_grid = h_grid) |>
    apply(MARGIN = 1, function(x) pracma::trapz(t_grid, x))


  Wm_max <- NW_max(curves = curves,
                   x = t_grid,
                   bw = h_grid)


  WN <- curves_select(curves = curves,
                      x = t_grid,
                      h = h_grid) |>
    aperm(c(2, 1, 3)) |>
    colSums()

  Nmu <- aperm(
    curves_select(curves = curves, x = t_grid, h = h_grid) * Wm_max,
    c(2, 1, 3)
  ) |>
    (\(x) (colSums(x) / WN**2)**(-1))()

  var_mu <- apply(param_smooth$sigma**2 / Nmu,
                  MARGIN = 2,
                  function(x) pracma::trapz(x = t_grid, y = x))

  reg_mu <- apply(param_smooth$moments2 * (1/WN - 1/length(curves)),
                  MARGIN = 2,
                  function(x) pracma::trapz(x = t_grid, y = x))

  risk_mu <- bias_mu + var_mu + reg_mu

  min_idx <- which.min(risk_mu)

  h_grid[min_idx]


}


#' Adaptive estimation of the mean function
#'
#' Estimates the mean function given a list of functional data, by adapting to
#' the regularity of the sample paths. Adaptation is obtained through bandwidth
#' selection, which is taken to minimise the mean integrated squared error (MISE)
#' of the mean function.
#'
#' @param Y_list List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param xout Vector containing the evaluation points for smoothing curves.
#' @param param_list List, containing the following estimated parameters, for
#' example from `estimate_parameters_FPCA`.
#' @param bw_grid Vector containing the grid of bandwidths to optimise over.
#' @returns Vector containing the estimated mean function.
#' @export

mu_ISE <- function(Y_list, xout, param_list, bw_grid) {

  if(length(param_list$t) != length(xout)) {
    stop("param_list must be on the same grid as grid_smooth!")
  }

  bw_mu <- bw_mean(curves = Y_list,
                   param_smooth = param_list,
                   h_grid = bw_grid,
                   t_grid = xout)

  X_hat <- smooth_curves(data = Y_list,
                         grid = xout,
                         bandwidth = bw_mu)$smoothed_curves |>
    abind::abind(along = 2)

  wi <- curves_select(curves = Y_list, x = xout, h = bw_mu)

  rowSums(x = wi * X_hat) / rowSums(wi)

  # list(
  #   t = xout,
  #   mu = rowSums(x = wi * X_hat) / rowSums(wi)
  # )

}




