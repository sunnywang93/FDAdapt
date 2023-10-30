
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
#' @param para_smooth List of estimated parameters.
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
#' @references Wang S., Patilea V., Klutchnikoff N., (2023+) - Adaptive
#' Functional Principal Components Analysis
#' @export


bw_FPCA <- function(data, param_smooth, h_grid, t_grid, psi_mat,
                    lambda_vec, inflate = TRUE) {


  # Compute constant kernel
  cst_kernel <- 1.5 * (1 / (1 + 2*param_smooth$H) - 1 / (3 + 2*param_smooth$H))

  # Compute bias term of eigenvalues ==========================================
  bias_rate_val <- bias_kernel(H = param_smooth$H,
                               L = param_smooth$L,
                               cst = cst_kernel,
                               h_grid = h_grid)


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
  bias_ss <- sapply(seq_len(ncol(psi_mat)), function(j) {
    sapply(seq_along(h_grid), function(h) {
      apply(bias_rate_val[h, ] * psi_mat[, -j]**2,
            2,
            function(x) pracma::trapz(t_grid, x))
    }, simplify = "array")
  }, simplify = "array")


  bias_st <- apply(psi_mat**2, 2, function(x) param_smooth$moments2 * x) |>
    apply(2, function(x) pracma::trapz(t_grid, x))

  bias_s_fun <- sweep(bias_ss,
                      3,
                      bias_st,
                      FUN = "*")

  bias_tt <- apply(psi_mat**2,
                   2,
                   function(x)
                     apply(x * bias_rate_val,
                           1,
                           function(y) pracma::trapz(t_grid, y)
                     )
                   )

  bias_ts <- sapply(seq_len(ncol(psi_mat)), function(j) {
    apply(psi_mat[, -j]**2 * param_smooth$moments2,
          2,
          function(x) pracma::trapz(t_grid, x)
    )
  }, simplify = "array")

  bias_t_fun <- sapply(seq_along(h_grid), function(h) {
    sweep(bias_ts,
          2,
          bias_tt[h, ],
          FUN = "*"
          )
  }, simplify = "array") |>
    aperm(c(1, 3, 2))


  # Compute bias constants
  fun_cst <- sapply(seq_along(lambda_vec),
                     function(j) (lambda_vec[j] - lambda_vec[-j])**(-2))

  bias_term_fun <- sapply(seq_len(ncol(psi_mat)), function(j) {
    colSums(fun_cst[, j] * (bias_t_fun[,, j] + bias_s_fun[,, j]))
  }, simplify = "array")

  # End of bias computation for eigenfunctions ================================

  # Compute variance term for eigenvalues =====================================
  # Compute variance at s
  var_s_num_val <- tcrossprod(param_smooth$sigma**2, param_smooth$moments2)

  # Compute Ngamma
  Ngamma <- N_gamma(curves = data,
                    x = t_grid,
                    bw = h_grid)


  # Compute variance term
  var_term_val <- sapply(seq_len(ncol(psi_mat)), function(j) {
    sapply(seq_along(h_grid), function(h) {
      (var_s_num_val / Ngamma$Ngamma[,, h] +
         t(var_s_num_val) / t(Ngamma$Ngamma[,,h])) *
        tcrossprod(psi_mat[, j]**2)
    }, simplify = "array")
  }, simplify = "array") |>
    apply(c(2, 3, 4), function(x) pracma::trapz(t_grid, x)) |>
    apply(c(2, 3), function(x) pracma::trapz(t_grid, x))

  # End of computation for eigenvalue variance =================================

  # Compute variance for eigenfunctions =======================================
  var_tst <- psi_mat**2 * param_smooth$moments2

  var_tss <- sapply(seq_len(ncol(psi_mat)), function(j) {
    psi_mat[, -j]**2 * param_smooth$sigma**2
  }, simplify = "array")

  var_ts <- sapply(seq_len(ncol(psi_mat)), function(j) {
    sapply(seq_len(ncol(psi_mat) - 1), function(k) {
      sapply(seq_along(h_grid), function(h) {
        outer(var_tss[, k, j], var_tst[, j]) / Ngamma$Ngamma[,, h]
      }, simplify = "array")
    }, simplify = "array")
  }, simplify = "array") |>
    apply(c(2, 3, 4, 5), function(x) pracma::trapz(t_grid, x)) |>
    apply(c(2, 3, 4), function(x) pracma::trapz(t_grid, x))

  var_stt <- psi_mat**2 * param_smooth$sigma**2

  var_sts <- sapply(seq_len(ncol(psi_mat)), function(j) {
    psi_mat[, -j]**2 * param_smooth$moments2
  }, simplify = "array")

  var_st <- sapply(seq_len(ncol(psi_mat)), function(j) {
    sapply(seq_len(ncol(psi_mat) - 1), function(k) {
      sapply(seq_along(h_grid), function(h) {
        outer(var_sts[, k, j], var_stt[, j]) / t(Ngamma$Ngamma[,, h])
      }, simplify = "array")
    }, simplify = "array")
  }, simplify = "array") |>
    apply(c(2, 3, 4, 5), function(x) pracma::trapz(t_grid, x)) |>
    apply(c(2, 3, 4), function(x) pracma::trapz(t_grid, x))

  var_term_fun <- sapply(seq_along(h_grid), function(h) {
    colSums((var_ts[h,,] + var_st[h,,]) * fun_cst)
  }) |>
    t()
  # End of variance term calculations for eigenfunctions ======================

  # Compute regularising term for eigenvalues
  reg_term_val <- sapply(seq_len(ncol(psi_mat)), function(j) {
    sapply(seq_along(h_grid), function(h) {
      param_smooth$moments2prod *
        ((1 / Ngamma$WN_bi[,, h]) - (1 / length(data))) *
        tcrossprod(psi_mat[, j]**2)
    }, simplify = "array")
  }, simplify = "array") |>
    apply(c(2, 3, 4), function(x) pracma::trapz(t_grid, x)) |>
    apply(c(2, 3), function(x) pracma::trapz(t_grid, x))

  # Compute regularising term for eigenfunctions
  reg_rate_fun <- sapply(seq_len(ncol(psi_mat)), function(j) {
    sapply(seq_len(ncol(psi_mat) - 1), function(k) {
      sapply(seq_along(h_grid), function(h) {
        outer(psi_mat[, j]**2, psi_mat[, -j][, k]**2) *
          param_smooth$moments2prod *
          ((1 / Ngamma$WN_bi[,, h]) - (1 / length(data)))
      }, simplify = "array")
    }, simplify = "array")
  }, simplify = "array") |>
    apply(c(2, 3, 4, 5), function(x) pracma::trapz(t_grid, x)) |>
    apply(c(2, 3, 4), function(x) pracma::trapz(t_grid, x))

  reg_term_fun <- sapply(seq_along(h_grid), function(h) {
    colSums(reg_rate_fun[h, ,] * fun_cst)
  }) |>
    t()

  # Compute risk
  risk_val <- 4 * bias_term_val + 2 * var_term_val + reg_term_val
  risk_fun <- 2 * bias_term_fun + 2 * var_term_fun + reg_term_fun

  min_idx_val <- apply(risk_val, 2, which.min)
  min_idx_fun <- apply(risk_fun, 2, which.min)

  # Obtain h*
  h_star_val <- h_grid[min_idx_val]
  h_star_fun <- h_grid[min_idx_fun]

  # Inflate by discretisation error if desired
  if(inflate) {
    h_constant_val <- log(1 / h_star_val)
    h_constant_fun <- log(1 / h_star_fun)
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
#' @param sigma Vector, containing the estimated noise for the purposes of
#' diagonal correction.
#' @param diag_cor Boolean, indicating whether to perform diagonal correction.
#' @returns Matrix of dimension `length(t_smooth) x length(t_smooth)`.
#' @references Wang S., Patilea V., Klutchnikoff N., (2023+) - Adaptive
#' Functional Principal Components Analysis
#' @export

cov_FPCA <- function(data, h, h_power, t_smooth, center = TRUE, sigma,
                     diag_cor) {

# Smooth curves at optimal bandwidth
smoothed <- smooth_curves(data = data, grid = t_smooth, bandwidth = h)

# Compute wi(t)
wi <- curves_select(curves = data,
                    x = t_smooth,
                    h = h)


WN <- rowSums(wi)

WN_bi <- apply(wi,
               2,
               tcrossprod,
               simplify = FALSE) |>
  (\(x) Reduce('+', x))()


# Center curves
if(center) {
  mu <- purrr::imap(smoothed$smoothed_curves,
                    ~.x * wi[, .y]) |>
    (\(x) Reduce('+', x) / WN)()

  smoothed_curves <- sapply(smoothed$smoothed_curves,
                                     function(x) x - mu)

  Gamma <- apply(wi * smoothed_curves,
                 2,
                 tcrossprod,
                 simplify = FALSE) |>
    (\(x) Reduce('+', x) / WN_bi)()

} else {
  Gamma <- purrr::imap(smoothed$smoothed_curves,
                       ~tcrossprod(wi[, .y] * .x)) |>
    (\(x) Reduce('+', x) / WN_bi)()
}


# Compute diagonal band correction

  if(diag_cor) {
    # Compute multiplicative term
    norm_diag <- tcrossprod(sigma) / WN_bi

    # Compute diagonal sum
    diag_sum <- purrr::imap(smoothed$weights,
                            ~crossprod(.x) * wi[, .y]) |>
      (\(x) Reduce('+', x))()

    # Correct diagonal band
    Gamma - norm_diag * diag_sum
  } else {
    Gamma
  }

}

#' Performs functional principal components analysis
#'
#' Perform FPCA of functional data adaptively using the local regularity of
#' curves. See references for more details on the algorithm.
#'
#' @param data List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param grid_smooth Vector of sampling points to smooth curves. Eigenfunctions
#' will be computed at these sampling points.
#' @param grid_bw Vector of points for bandwidth optimization.
#' @param param_list List, containing the estimated parameters.
#' @param inflate_bw Boolean, indicating whether to correct the bandwidth due
#' to discretization error from estimating the regularity.
#' @param center Boolean, indicating whether to center curves.
#' @param nelements Numeric, indicating the number of eigen-elements to keep.
#' @param h_power Numeric, power to raise the bandwidth to when estimating
#' sigma for the purposes of diagonal correction.
#' @returns List, containing the following elements:
#' - **$params** List containing the estimated parameters.
#' - **$bw** Numeric containing the bandwidth used for smoothing curves.
#' - **$evalues** Vector containing the normalised eigenvalues.
#' - **$efunctions** Matrix containing the normalised eigenfunctions, with the j-th
#' column representing the j-th eigenfunction.
#' @export

FPCA <- function(data,
                 grid_smooth,
                 grid_bw,
                 param_list,
                 inflate_bw = TRUE,
                 center = TRUE,
                 nelements = 10,
                 diag_cor = TRUE,
                 h_power = 0.9) {


  # Compute covariance function
  covariance_FPCA <- cov_ISE(Y_list = data,
                             xout = grid_smooth,
                             hout = grid_bw,
                             param_list = param_list,
                             h_power = h_power,
                             inflate_bw = TRUE,
                             diag_cor = diag_cor)


  # Perform eigen-analysis and return normalized elements
  elements <- normalise_eigen(covariance_FPCA$gamma,
                              nelements = nelements)

  list(params = param_list,
       bw = covariance_FPCA$bw,
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
#' @param param_list List of estimated parameters, for example from the
#' function `estimate_parameters_FPCA`. Must be on the same grid as `grid_smooth`,
#' which can be interpolated for example using `intp_list`.
#' @param inflate_bw Boolean, indicating whether to correct the bandwidth due
#' to discretization error from estimating the regularity.
#' @param center Boolean, indicating whether to center curves.
#' @param nelements Numeric, indicating the number of eigen-elements to keep.
#' Numeric must be at least 3.
#' @param diag_cor Boolean, indicating whether to perform diagonal correction.
#' @param h_power Numeric, power to raise bandwidth to when estimating sigma
#' for the purposes of diagonal correction for the covariance function.
#' @returns List, containing the following elements:
#' - **$evalues** Vector containing the normalised eigenvalues.
#' - **$efunctions** Matrix containing the normalised eigenfunctions, with the j-th
#' column representing the j-th eigenfunction.
#' - **$bw_val** Vector containing the bandwidth used to smooth curves for each
#' eigenvalue.
#' - **$bw_fun** Vector containing the bandwidth used to smooth curves for each
#' eigenfunction.
#' @export


FPCA_adapt <- function(data, grid_smooth, grid_bw, param_list,
                       inflate_bw = TRUE, center = TRUE, nelements = 10,
                       diag_cor = TRUE, h_power = 0.9) {

  if(length(param_list$t) != length(grid_smooth)) {
    stop("param_list must be on the same grid as grid_smooth!")
  }

  # Obtain preliminary estimates of eigen-elements
  eelements_prelim <- FPCA(data = data,
                           grid_smooth = grid_smooth,
                           grid_bw = grid_bw,
                           param_list = param_list,
                           inflate_bw = inflate_bw,
                           center = center,
                           nelements = nelements,
                           diag_cor = diag_cor,
                           h_power = h_power)

  # Obtain adaptive bandwidth estimates
  bw_adaptive <- bw_FPCA(data = data,
                         param_smooth = param_list,
                         h_grid = grid_bw,
                         t_grid = grid_smooth,
                         psi_mat = eelements_prelim$efunctions,
                         lambda_vec = eelements_prelim$evalues,
                         inflate = inflate_bw)

  # Get covariance of relevant bandwidths and eigen-elements
  eigen_val <- purrr::map(bw_adaptive$bw_val, ~cov_FPCA(data = data,
                                                        h = .x,
                                                        h_power = h_power,
                                                        t_smooth = grid_smooth,
                                                        center = center,
                                                        sigma = param_list$sigma,
                                                        diag_cor = diag_cor)) |>
    purrr::map(~normalise_eigen(.x, nelements = nelements)$values)

  val_adapt <- purrr::map_dbl(seq_len(length(eigen_val)), ~eigen_val[[.x]][.x])

  eigen_fun <- purrr::map(bw_adaptive$bw_fun, ~cov_FPCA(data = data,
                                                        h = .x,
                                                        h_power = h_power,
                                                        t_smooth = grid_smooth,
                                                        center = center,
                                                        sigma = param_list$sigma,
                                                        diag_cor = diag_cor)) |>
    purrr::map(~normalise_eigen(.x, nelements = nelements)$vectors)

  fun_adapt <- sapply(seq_len(length(eigen_fun)),
                      function(j) eigen_fun[[j]][, j])

  # Apply Gram-Schmidt to ensure orthonormality of eigenfunctions
  fun_ortho <- ortho_funmat(X = fun_adapt,
                            t = grid_smooth)

  list(evalues = val_adapt,
       efunctions = fun_ortho,
       bw_val = bw_adaptive$bw_val,
       bw_fun = bw_adaptive$bw_fun)

}





