#' Computes PACE estimates with learning and online set
#'
#' Adaptive PACE estimates are computed given a list of learning curves and
#' a list of online curves. Quantities are estimated on the set of learning
#' curves, before being interpolated onto the sampling points given by a list
#' of online curves. Each quantity is adaptively and separately computed
#' using optimal estimates in the sense of mean-integrated squared error (MISE).
#'
#' @param Y_learn List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param Y_online List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param x_out Vector containing the evaluation points to smooth curves.
#' @param h_out Vector containing the grid of bandwidths.
#' @param param_out Vector containing the evaluation points for parameter
#' estimation.
#' @param ncv_param Numeric containing the number of curves used for learning
#' the cross-validation bandwidth in presmoothing.
#' @param h_q Numeric containing the quantile selected from the set of
#' cross-validation bandwidths.
#' @param gamma_H Numeric, indicating the gamma to be used for estimation of H.
#' See `estimate_regularity`.
#' @param gamma_L Numeric, indicating the gamma to be used for estimation of L.
#' See `estimate_regularity`.
#' @param n_knots Numeric indicating the number of knots used for smoothing H and L,
#' where splines are used.
#' @param h_power Numeric indicating the power to raise bandwidth to when
#' estimating sigma for the purposes of diagonal correction for the covariance function.
#' @param n_scores Numeric indicating the number of scores to estimate.
#' @returns Matrix of dimension K x N, where K is the number of scores and
#' N is the number of curves in `Y_online`.
#' @export



PACE_adapt <- function(Y_learn, Y_online, x_out, h_out, param_list,
                       h_power, n_scores
                       ) {


  elements_learn <- FPCA_adapt(data = Y_learn,
                               grid_smooth = x_out,
                               grid_bw = h_out,
                               h_power = h_power,
                               true_params = param_list)


  mu_learn <- mu_ISE(Y_list = Y_learn,
                     xout = x_out,
                     param_list = param_list,
                     bw_grid = h_out)


  cov_learn <- cov_ISE(Y_list = Y_learn,
                       xout = x_out,
                       hout = h_out,
                       param_list = param_list,
                       h_power = h_power,
                       inflate_bw = TRUE)


  v_online <- purrr::map(Y_online,
                         ~apply(elements_learn$efunctions,
                                2,
                                function(psi) interpolate1D(x = x_out,
                                                            y = psi,
                                                            xout = .x$t)$y
                                ))


  mu_online <- purrr::map(Y_online,
                          ~interpolate1D(x = x_out,
                                         y = mu_learn,
                                         xout = .x$t)$y)


  noise <- estimate_sigma(data = Y_online,
                          sigma_grid = xout)

  noise_online <- purrr::map(Y_online,
                             ~interpolate1D(x = x_out,
                                            y = noise,
                                            xout = .x$t)$y)

  cov_online <- purrr::map2(Y_online, noise_online,
                           ~interpolate2D(x = x_out,
                                          y = x_out,
                                          z = cov_learn,
                                          xout = .x$t,
                                          yout = .x$t) + diag(.y))

  Y_cen <- purrr::map2(Y_online, mu_online,
                       ~.x$x - .y)

  v_prod <- purrr::map(v_online,
                       ~sweep(.x, 2, elements_learn$evalues, FUN = "*"))


  sapply(seq_along(Y_online), function(i) {
    apply(v_prod[[i]], 2,
          function(vk) t(vk) %*% solve(cov_online[[i]]) %*% Y_cen[[i]])
  })


}









