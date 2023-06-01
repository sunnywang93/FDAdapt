#' Compute FPCA bandwidth without index specific constants
#'
#' Performs estimation of the FPCA bandwidth when all the *index specific*
#' constants are dropped for the risk bound, resulting in a single bandwidth
#' for both the eigenvalues and the eigenfunctions, instead of a different
#' bandwidth for each eigenvalue and eigenfunction.
#'
#' @param data List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param parameters List of parameters, which can be obtained from
#' `estimate_parameters_FPCA`.
#' @param h_grid Vector of points for bandwidth optimization.
#' @param t_grid Vector of sampling points to smooth curves.
#' @param inflate Boolean, indicating whether to inflate bandwidth as a
#' correction for the discretization error from estimating regularity.
#' @returns Numeric, containing the optimal bandwidth.
#' @references Wang S., Patilea V., Klutchnikoff N., (2023+) - Adaptive
#' Functional Principal Components Analysis
#' @export


bw_FPCA <- function(data, parameters, h_grid, t_grid, inflate = TRUE,
                    interp_type = "linear") {

  # Interpolate 1D-parameters onto smoothing grid
  if(length(parameters$t) != length(t_grid)) {
    param_smooth <- lapply(parameters[c("H", "L", "sigma", "moments2")],
                           function(i) interpolate1D(x = parameters$t,
                                                     y = i,
                                                     xout = t_grid,
                                                     type = interp_type)$y)

  # Interpolate 2D-parameter
  mom2prod_smooth <- interpolate2D(x = parameters$t,
                                   y = parameters$t,
                                   z = parameters$moments2prod,
                                   xout = t_grid,
                                   yout = t_grid)

  # Append to list of parameters
  param_smooth$moments2prod <- mom2prod_smooth

  } else {
    param_smooth <- parameters
  }

  # Compute constant kernel
  cst_kernel <- 1.5 * (1 / (1 + 2*param_smooth$H) - 1 / (3 + 2*param_smooth$H))

  # Compute bias rate term
  bias_t <- outer(h_grid, 2 * param_smooth$H, FUN = "^") |>
    sweep(MARGIN = 2, STATS = param_smooth$L**2, FUN = "*") |>
    sweep(MARGIN = 2, STATS = cst_kernel, FUN = "*") |>
    apply(MARGIN = 1, FUN = function(x) pracma::trapz(t_grid, x))

  # Compute integrated bias at s
  bias_s <- pracma::trapz(t_grid, param_smooth$moments2)

  # Compute bias term
  bias_term <-  4 * bias_t * bias_s

  # Compute variance at (s|t)
  variance_s_num <- outer(param_smooth$sigma**2, param_smooth$moments2)

  # Compute sample weights
  wi <- purrr::map(data, ~abs(outer(.x$t, t_grid, FUN = "-"))) |>
    purrr::map(~outer(.x, h_grid, FUN = "<=") * 1) |>
    purrr::map(~(colSums(.x) >= 1) * 1) |>
    abind::abind(along = 3) |>
    aperm(c(1, 3, 2))

  WN <- apply(wi, c(1, 3), sum)

  WN_bi <- apply(wi, c(2, 3), function(x) tcrossprod(x)) |>
    array(dim = c(length(t_grid),
                  length(t_grid),
                  length(data),
                  length(h_grid))) |>
    apply(c(1, 2, 4), sum)

  # Compute kernel weights
  weights_max <- lapply(data, function(i) {
    epa_kernel(outer(outer(i$t, t_grid, FUN = "-"),
                     h_grid, FUN = "/")) |>
      apply(c(2, 3), max)
  })

  weights_denom <- lapply(data, function(i) {
    epa_kernel(outer(outer(i$t, t_grid, FUN = "-"),
                     h_grid, FUN = "/")) |>
      colSums()
  })

  weights_max_norm <- purrr::map2(weights_max, weights_denom, ~(.x / .y)) |>
    rapply(f = function(x) ifelse(is.nan(x), 0, x), how = "replace") |>
    abind::abind(along = 3) |>
    aperm(c(1, 3, 2))

  # Compute Ngamma's
  Ngamma <- sapply(seq_along(h_grid), function(h) {
    sapply(seq_along(data), function(i) {
      sweep(tcrossprod(wi[, i, h]), 2, weights_max_norm[, i, h], FUN = "*")
    }, simplify = "array")
  }, simplify = "array") |>
    apply(c(1,2,4), sum) |>
    (\(x) 1 / (x / WN_bi**2))()

  # Compute variance term
  variance_term <- sapply(seq_along(h_grid), function(h) {
    variance_s <- variance_s_num / Ngamma[,, h]
    variance_t <- t(variance_s_num) / t(Ngamma[,, h])
    variance_s + variance_t
  }, simplify = "array") |>
    apply(c(2, 3), function(x) pracma::trapz(t_grid, x)) |>
    apply(2, function(x) pracma::trapz(t_grid, x))

  # Compute regularising term
  regularising_term <- sapply(seq_along(h_grid), function(h) {
    ((1 / WN_bi[,,h]) - (1 / length(data))) * param_smooth$moments2prod
  }, simplify = "array") |>
    apply(c(2, 3), function(x) pracma::trapz(t_grid, x)) |>
    apply(2, function(x) pracma::trapz(t_grid, x))

  # Compute risk
  risk <- bias_term + 2 * variance_term + regularising_term

  min_idx <- which.min(risk)

  # Obtain h*
  h_star <- h_grid[min_idx]

  # Inflate by discretisation error if desired
  if(inflate) {
    m <- mean(purrr::map_dbl(data, ~length(.x$t)))
    h_constant <- log(1 / h_star)**(log(1 / h_star) / log(length(data) * m))
    h_star <- h_constant * h_star
  }
  # Return h_star
  h_star

}

#' Compute fully adaptive FPCA bandwidth
#'
#' Performs estimation of the FPCA bandwidth for each index of the eigenvalues
#' and eigenfunctions separately, leading to 2 * J different bandwidths, where
#' J indicates the number of eigen-elements considered. Differences in bandwidths
#' are driven by the different constants specific to each eigen-element.
#'
#' @param data List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param parameters List of parameters, which can be obtained from
#' `estimate_parameters_FPCA`.
#' @param h_grid Vector of points for bandwidth optimization.
#' @param t_grid Vector of sampling points to smooth curves, should be in the
#' interval `[0, 1]`.
#' @param psi_mat Matrix containing the eigenfunctions, with columns
#' representing the j-th eigenfunction and the rows are evaluated points on
#' a grid.
#' @param lambda_vec Vector of eigenvalues, whose length must match the number
#' of columns of `psi_mat`.
#' @param inflate Boolean, indicating whether to inflate bandwidth as a
#' correction for the discretization error from estimating regularity.
#' @returns Numeric, containing the optimal bandwidth.
#' @references Wang S., Patilea V., Klutchnikoff N., (2023+) - Adaptive
#' Functional Principal Components Analysis
#' @export


bw_FPCA_adapt <- function(data, parameters, h_grid, t_grid, psi_mat,
                          lambda_vec, inflate = TRUE, interp_type = "linear") {

  # Interpolate 1D-parameters onto smoothing grid
  if(length(parameters$t) != length(t_grid)) {
    param_smooth <- lapply(parameters[c("H", "L", "sigma", "moments2")],
                           function(i) interpolate1D(x = parameters$t,
                                                     y = i,
                                                     xout = t_grid,
                                                     type = interp_type)$y)

  # Interpolate 2D-parameter
    mom2prod_smooth <- interpolate2D(x = parameters$t,
                                     y = parameters$t,
                                     z = parameters$moments2prod,
                                     xout = t_grid,
                                     yout = t_grid)

  # Append to list of parameters
    param_smooth$moments2prod <- mom2prod_smooth
  } else {
    param_smooth <- parameters
  }

  # Compute constant kernel
  cst_kernel <- 1.5 * (1 / (1 + 2*param_smooth$H) - 1 / (3 + 2*param_smooth$H))

  # Compute bias term of eigenvalues ==========================================
  bias_rate_val <- outer(h_grid, 2 * param_smooth$H, FUN = "^") |>
    sweep(MARGIN = 2, STATS = param_smooth$L**2, FUN = "*") |>
    sweep(MARGIN = 2, STATS = cst_kernel,  FUN = "*")

  bias_t_val <- apply(psi_mat**2, 2,
                      function(x) sweep(bias_rate_val, 2, x, FUN = "*") |>
                        apply(1, function(y) pracma::trapz(t_grid, y)))

  # Compute integrated bias at s
  bias_s_val <- apply(psi_mat**2, 2,
                      function(x) pracma::trapz(t_grid,
                                                param_smooth$moments2 * x))

  # Compute bias term
  bias_term_val <- sweep(bias_t_val, 2, bias_s_val, FUN = "*")

  #End of bias term for eigenvalues ===========================================

  # Compute bias terms of eigenfunctions ======================================
  bias_rate_ss <- outer(h_grid, 2 * param_smooth$H, FUN = "^") |>
    sweep(MARGIN = 2, STATS = param_smooth$L**2, FUN = "*") |>
    sweep(MARGIN = 2, STATS = cst_kernel,  FUN = "*")

  bias_ss <- sapply(seq_len(ncol(psi_mat)), function(j) {
    sweep(bias_rate_ss, 2, rowSums(psi_mat[, -j]**2), FUN = "*") |>
      apply(1, function(s) pracma::trapz(t_grid, s))
  }, simplify = "array")

  bias_st <- apply(psi_mat, 2, function(x) param_smooth$moments2 * x**2) |>
    apply(2, function(x) pracma::trapz(t_grid, x))

  bias_s_fun <- sweep(bias_ss, 2, bias_st, FUN = "*")

  bias_tt <- sapply(seq_len(ncol(psi_mat)), function(j) {
    sweep(bias_rate_ss, 2, psi_mat[, j]**2, FUN = "*") |>
      apply(1, function(t) pracma::trapz(t_grid, t))
  }, simplify = "array")

  bias_ts <- sapply(seq_len(ncol(psi_mat)), function(j) {
    pracma::trapz(t_grid, rowSums(psi_mat[, -j]**2) * param_smooth$moments2)
  })

  bias_t_fun <- sweep(bias_tt, 2, bias_ts, FUN = "*")

  # Compute bias constants
  bias_cst <- sapply(seq_along(lambda_vec),
                     function(j) 1 / sum((lambda_vec[j] - lambda_vec[-j])**2))

  bias_term_fun <- sweep(bias_t_fun + bias_s_fun, 2, bias_cst, FUN = "*")

  # End of bias computation for eigenfunctions ================================

  # Compute variance term for eigenvalues =====================================
  # Compute variance at s
  variance_s_num_val <- tcrossprod(param_smooth$sigma**2, param_smooth$moments2)

  # Compute sample weights
  wi <- purrr::map(data, ~abs(outer(.x$t, t_grid, FUN = "-"))) |>
    purrr::map(~outer(.x, h_grid, FUN = "<=") * 1) |>
    purrr::map(~(colSums(.x) >= 1) * 1)

  WN <- Reduce('+', wi)

  WN_bi <- sapply(seq_along(h_grid), function(h) {
    Reduce('+',
      lapply(wi, function(i) tcrossprod(i[, h]))
    )
  }, simplify = "array")

  # Compute kernel weights
  weights_max <- purrr::map(data, ~outer(.x$t, t_grid, "-")) |>
    purrr::map(~apply(epa_kernel(outer(.x, h_grid, FUN = "/")), c(2, 3), max))

  weights_denom <- purrr::map(data, ~outer(.x$t, t_grid, FUN = "-")) |>
    purrr::map(~colSums(epa_kernel(outer(.x, h_grid, FUN = "/"))))

  weights_max_norm <- purrr::map2(weights_max, weights_denom, ~.x / .y) |>
    rapply(f = function(x) ifelse(is.nan(x), 0, x), how = "replace")

  # Compute Ngamma's (unnormalised)
  Ngamma <- sapply(seq_along(h_grid), function(h) {
    Reduce('+',
      lapply(seq_along(data), function(i) {
        tcrossprod(wi[[i]][, h]) * weights_max_norm[[i]][, h]
     })
    )
  }, simplify = "array")

  Ngamma <- 1 / (Ngamma / WN_bi**2)


  # Compute variance term
  variance_term_val <- sapply(seq_len(ncol(psi_mat)), function(j) {
    sapply(seq_along(h_grid), function(h) {
      variance_s <- (variance_s_num_val * tcrossprod(psi_mat[, j]**2)) / Ngamma[,,h]
      variance_t <- (t(variance_s_num_val) * t(tcrossprod(psi_mat[,j]**2))) /
        t(Ngamma[,, h])
      variance_s + variance_t
    }, simplify = "array") |>
    apply(c(2, 3), function(x) pracma::trapz(t_grid, x)) |>
      apply(2, function(x) pracma::trapz(t_grid, x))
  })

  # End of computation for eigenvalue variance =================================

  # Compute variance for eigenfunctions =======================================
  variance_st <- apply(psi_mat**2, 2, function(x) x * param_smooth$moments2)

  variance_ss <- sapply(seq_len(ncol(psi_mat)), function(j) {
    rowSums(psi_mat[, -j]**2) * param_smooth$sigma**2
  })

  variance_s_fun <- sapply(seq(ncol(psi_mat)), function(j) {
    sapply(seq_along(h_grid), function(h) {
      tcrossprod(variance_st[, j], variance_ss[, j]) / Ngamma[,, h]
    }, simplify = "array")
  }, simplify = "array") |>
    apply(c(2, 3, 4), function(x) pracma::trapz(t_grid, x)) |>
    apply(c(2, 3), function(x) pracma::trapz(t_grid, x))

  variance_tt <- apply(psi_mat**2, 2, function(x) x * param_smooth$sigma**2)

  variance_ts <- sapply(seq_len(ncol(psi_mat)), function(j) {
    rowSums(psi_mat[, -j]**2) * param_smooth$moments2
  })

  variance_t_fun <- sapply(seq(ncol(psi_mat)), function(j) {
    sapply(seq_along(h_grid), function(h) {
      tcrossprod(variance_tt[, j], variance_ts[, j]) / t(Ngamma[,, h])
    }, simplify = "array")
  }, simplify = "array") |>
    apply(c(2, 3, 4), function(x) pracma::trapz(t_grid, x)) |>
    apply(c(2, 3), function(x) pracma::trapz(t_grid, x))

  variance_cst <- sapply(seq_along(lambda_vec),
                     function(j) 1 / sum((lambda_vec[j] - lambda_vec[-j])**2))


  variance_term_fun <- sweep(variance_t_fun + variance_s_fun,
                             2, variance_cst, FUN = "*")

  # End of variance term calculations for eigenfunctions ======================

  # Compute regularising term for eigenvalues
  regularising_term_val <- sapply(seq_len(ncol(psi_mat)), function(j) {
    sapply(seq_along(h_grid), function(h) {
      ((1 / WN_bi[,, h]) - (1 / length(data))) * param_smooth$moments2prod *
        tcrossprod(psi_mat[, j]**2)
    }, simplify = "array") |>
      apply(c(2, 3), function(x) pracma::trapz(t_grid, x)) |>
      apply(2, function(x) pracma::trapz(t_grid, x))
  })

  # Compute regularising term for eigenfunctions
  regularising_term_fun <- sapply(seq_len(ncol(psi_mat)), function(j) {
    sapply(seq_along(h_grid), function(h) {
      tcrossprod(psi_mat[, j]**2, rowSums(psi_mat[, -j]**2)) *
        param_smooth$moments2prod *
        ((1 / WN_bi[,, h]) - (1 / length(data)))
    }, simplify = "array") |>
      apply(c(2, 3), function(x) pracma::trapz(t_grid, x)) |>
      apply(2, function(x) pracma::trapz(t_grid, x))
  }) |>
    sweep(2, variance_cst, FUN = "*")


  # Compute risk
  risk_val <- 4 * bias_term_val + 2 * variance_term_val + regularising_term_val
  risk_fun <- 2 * bias_term_fun + 2 * variance_term_fun + regularising_term_fun

  min_idx_val <- apply(risk_val, 2, which.min)
  min_idx_fun <- apply(risk_fun, 2, which.min)

  # Obtain h*
  h_star_val <- h_grid[min_idx_val]
  h_star_fun <- h_grid[min_idx_fun]

  # Inflate by discretisation error if desired
  if(inflate) {
    m <- mean(purrr::map_dbl(data, ~length(.x$t)))
    h_constant_val <- log(1 / h_star_val)**(log(1 / h_star_val) /
                                              log(length(data) * m))
    h_constant_fun <- log(1 / h_star_fun)**(log(1 / h_star_fun) /
                                              log(length(data) * m))
    h_star_val <- h_constant_val * h_star_val
    h_star_fun <- h_constant_fun * h_star_fun
  }
  # Return h_star
  list(bw_val= h_star_val,
       bw_fun = h_star_fun)

}




#' Compute the FPCA covariance function
#'
#' Estimates the covariance function for the purposes of FPCA.
#'
#' @param data List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param h Numeric containing the optimal bandwidth.
#' @param h_power Numeric, power to raise the bandwidth to when estimating
#' sigma for the purposes of diagonal correction.
#' @param t_smooth Vector of sampling points to smooth curves. Covariance will
#' be computed at bi-dimensional grid points constructed from this vector.
#' @param center Boolean, indicating whether to center curves by their mean
#' function.
#' @returns Matrix of dimension `length(t_smooth) x length(t_smooth)`.
#' @references Wang S., Patilea V., Klutchnikoff N., (2023+) - Adaptive
#' Functional Principal Components Analysis
#' @export

cov_FPCA <- function(data, h, h_power, t_smooth, center = TRUE, sigma_true,
                     diag_cor) {

# Smooth curves at optimal bandwidth
smoothed <- smooth_curves(data = data, grid = t_smooth, bandwidth = h)

# Compute wi(t)
wi <- purrr::map(data, ~(abs(outer(.x$t, t_smooth, FUN = "-")) <= h) * 1) |>
  purrr::map(~(colSums(.x) >= 1) * 1)

# Compute wi(s)
WN <- Reduce('+', wi)

# Center curves
if(center) {
  mu <- Reduce('+', purrr::map2(wi, smoothed$smoothed_curves, ~.x * .y)) / WN

  smoothed$smoothed_curves <- purrr::map(smoothed$smoothed_curves, ~.x - mu)
}

# Compute cov
Gamma <- Reduce('+', purrr::map2(smoothed$smoothed_curves, wi,
                                 ~tcrossprod(.x) * tcrossprod(.y)))

WN_bi <- Reduce('+', purrr::map(wi, ~tcrossprod(.x)))

Gamma_cen <- Gamma / WN_bi

# Compute diagonal band correction

# Compute sigma's for the diagonal band
  if(missing(sigma_true)) {
    sigma_hat <- estimate_sigma(data = data, sigma_grid = t_smooth, h = h,
                                h_power = h_power)
  } else {
    sigma_hat <- sigma_true
  }

  if(diag_cor) {
    # Compute multiplicative term
    norm_diag <- tcrossprod(sigma_hat) / WN_bi
    # Compute diagonal sum
    diag_sum <- Reduce('+',
                       purrr::map2(wi, smoothed$weights,
                                   ~tcrossprod(.x) * crossprod(.y)))
    # Compute diagonal bias
    diag_bias <- diag_sum * norm_diag
    # Correct diagonal band
    Gamma_cen - diag_bias
  } else {
    Gamma_cen
  }

}

#' Performs functional principal components analysis
#'
#' Perform FPCA of functional data adaptively using the local regularity of
#' curves. See references for more details on the algorith.
#'
#' @param data List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param grid_smooth Vector of sampling points to smooth curves. Eigenfunctions
#' will be computed at these sampling points.
#' @param grid_bw Vector of points for bandwidth optimization.
#' @param grid_param Vector of sampling points to estimate parameters.
#' @param cv_set Numeric, containing the number of curves used for learning
#' the cross-validation bandwidth in presmoothing.
#' @param quantile Numeric, containing the quantile selected from the set of
#' cross-validation bandwidths.
#' @param inflate_bw Boolean, indicating whether to correct the bandwidth due
#' to discretization error from estimating the regularity.
#' @param center Boolean, indicating whether to center curves.
#' @param interp_type String, indicating the type of interpolation to perform
#' for the parameters from `grid_param` onto `grid_smooth`. Options include
#' c("linear", "constant", "nearest", "spline", "cubic").
#' @param nelements Numeric, indicating the number of eigen-elements to keep.
#' @param gamma_H Numeric, indicating the gamma to be used for estimation of
#' `H`. See `estimate_regularity`.
#' @param gamma_L Numeric, indicating the gamma to be used for estimation of
#' `L`. See `estimate_regularity`.
#' @param n_knots Numeric, number of knots used for smoothing H and L, where
#' splines are used.
#' @param h_power Numeric, power to raise the bandwidth to when estimating
#' sigma for the purposes of diagonal correction.
#' @param intp_param Boolean, where `TRUE` indicates using presmoothing + interpolation
#' to estimate regularity parameters, as compared to only presmoothing.
#' @returns List, containing the following elements:
#' - **$params** List containing the estimated parameters.
#' - **$bw** Numeric containing the bandwidth used for smoothing curves.
#' - **$evalues** Vector containing the eigenvalues.
#' - **$efunctions** Matrix containing the eigenfunctions, with the j-th
#' column representing the j-th eigenfunction.
#' @export

FPCA <- function(data, grid_smooth, grid_bw, grid_param, cv_set,
                 quantile, inflate_bw = TRUE,
                 center = TRUE, interp_method = "linear",
                 nelements = 10, gamma_H, gamma_L, n_knots, h_power,
                 intp_param, true_params, diag_cor = TRUE) {

  # Estimate parameters
  if(missing(true_params)) {
    params_FPCA <- estimate_parameters_FPCA(data = data,
                                            grid_points = grid_param,
                                            n_learn = cv_set,
                                            h_quantile = quantile,
                                            intp = intp_param,
                                            gamma_H = gamma_H,
                                            gamma_L = gamma_L,
                                            nknots = n_knots)
  } else {
    params_FPCA <- true_params
  }

  # Estimate optimal bandwidth
  bandwidth_FPCA <- bw_FPCA(data = data,
                            parameters = params_FPCA,
                            h_grid = grid_bw,
                            t_grid = grid_smooth,
                            inflate = inflate_bw,
                            interp_type = interp_method)

  # Compute covariance function
  if(missing(true_params)) {
    covariance_FPCA <- cov_FPCA(data = data,
                                h = bandwidth_FPCA,
                                h_power = h_power,
                                t_smooth = grid_smooth,
                                center = center,
                                diag_cor = diag_cor)
  } else {
    covariance_FPCA <- cov_FPCA(data = data,
                                h = bandwidth_FPCA,
                                t_smooth = grid_smooth,
                                center = center,
                                sigma_true = true_params$sigma,
                                diag_cor = diag_cor)
  }


  # Perform eigen-analysis and return normalized elements
  elements <- normalise_eigen(covariance_FPCA, nelements = nelements)

  list(params = params_FPCA,
       bw = bandwidth_FPCA,
       evalues = elements$values,
       efunctions = elements$vectors)

}

#' Performs adaptive functional principal components analysis
#'
#' Perform FPCA of functional data adaptively using the local regularity of
#' curves, and with bandwidths specific to each eigenvalue and eigenfunction.
#' Although more computationally demanding, it is more precise with regards to
#' the element-specific bandwidth.
#'
#' @param data List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param grid_smooth Vector of sampling points to smooth curves. Eigenfunctions
#' will be computed at these sampling points.
#' @param grid_bw Vector of points for bandwidth optimization.
#' @param grid_param Vector of sampling points to estimate parameters.
#' @param cv_set Numeric, containing the number of curves used for learning
#' the cross-validation bandwidth in presmoothing.
#' @param quantile Numeric, containing the quantile selected from the set of
#' cross-validation bandwidths.
#' @param inflate_bw Boolean, indicating whether to correct the bandwidth due
#' to discretization error from estimating the regularity.
#' @param center Boolean, indicating whether to center curves.
#' @param interp_type String, indicating the type of interpolation to perform
#' for the parametrs from `grid_param` onto `grid_smooth`. Options include
#' c("linear", "constant", "nearest", "spline", "cubic").
#' @param nelements Numeric, indicating the number of eigen-elements to keep.
#' @param gamma_H Numeric, indicating the gamma to be used for estimation of
#' `H`. See `estimate_regularity`.
#' @param gamma_L Numeric, indicating the gamma to be used for estimation of
#' `L`. See `estimate_regularity`.
#' @param n_knots Numeric, number of knots used for smoothing H and L, where
#' splines are used.
#' @param h_power Numeric, power to raise bandwidth to when estimating sigma
#' for the purposes of diagonal correction for the covariance function.
#' @param intp_param Boolean, where `TRUE` indicates using presmoothing + interpolation
#' to estimate regularity parameters, as compared to only presmoothing.
#' @returns List, containing the following elements:
#' - **$evalues** Vector containing the eigenvalues.
#' - **$efunctions** Matrix containing the eigenfunctions, with the j-th
#' column representing the j-th eigenfunction.
#' - **$bw_val** Vector containing the bandwidth used to smooth curves for each
#' eigenvalue.
#' - **$bw_fun** Vector containing the bandwidth used to smooth curves for each
#' eigenfunction.
#' @export


FPCA_adapt <- function(data, grid_smooth, grid_bw, grid_param, cv_set,
                       quantile, inflate_bw = TRUE,
                       center = TRUE, interp_method = "linear",
                       nelements = 10, gamma_H, gamma_L, n_knots,
                       h_power, intp_param, true_params, diag_cor = TRUE) {

  # Obtain preliminary estimates of eigen-elements
  if(missing(true_params)) {
    eelements_prelim <- FPCA(data = data,
                             grid_smooth = grid_smooth,
                             grid_bw = grid_bw,
                             grid_param = grid_param,
                             cv_set = cv_set,
                             quantile = quantile,
                             inflate_bw = inflate_bw,
                             center = center,
                             interp_method = interp_method,
                             nelements = nelements,
                             gamma_H = gamma_H,
                             gamma_L = gamma_L,
                             n_knots = n_knots,
                             h_power = h_power,
                             intp_param = intp_param,
                             diag_cor = diag_cor)
  } else {
    eelements_prelim <- FPCA(data = data,
                             grid_smooth = grid_smooth,
                             grid_bw = grid_bw,
                             grid_param = grid_param,
                             cv_set = cv_set,
                             quantile = quantile,
                             inflate_bw = inflate_bw,
                             center = center,
                             interp_method = interp_method,
                             nelements = nelements,
                             true_params = true_params,
                             diag_cor = diag_cor)
  }


  # Obtain adaptive bandwidth estimates
  bw_adaptive <- bw_FPCA_adapt(data = data,
                               parameters = eelements_prelim$params,
                               h_grid = grid_bw,
                               t_grid = grid_smooth,
                               psi_mat = eelements_prelim$efunctions,
                               lambda_vec = eelements_prelim$evalues,
                               inflate = inflate_bw,
                               interp_type = interp_method)

  # Get covariance of relevant bandwidths and eigen-elements
  eigen_val <- purrr::map(bw_adaptive$bw_val, ~cov_FPCA(data = data,
                                                        h = .x,
                                                        h_power = h_power,
                                                        t_smooth = grid_smooth,
                                                        center = center,
                                                        diag_cor = diag_cor)) |>
    purrr::map(~normalise_eigen(.x, nelements = nelements)$values)

  val_adapt <- purrr::map_dbl(seq_len(length(eigen_val)), ~eigen_val[[.x]][.x])

  eigen_fun <- purrr::map(bw_adaptive$bw_fun, ~cov_FPCA(data = data,
                                                        h = .x,
                                                        h_power = h_power,
                                                        t_smooth = grid_smooth,
                                                        center = center,
                                                        diag_cor = diag_cor)) |>
    purrr::map(~normalise_eigen(.x, nelements = nelements)$vectors)

  fun_adapt <- sapply(seq_len(length(eigen_fun)),
                      function(j) eigen_fun[[j]][, j])

  list(evalues = val_adapt,
       efunctions = fun_adapt,
       bw_val = bw_adaptive$bw_val,
       bw_fun = bw_adaptive$bw_fun)

}





