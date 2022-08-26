#import dependent scripts
source("/Users/swang/Dropbox/funeigen/covariance/cov_optimised.R")
source("/Users/swang/Dropbox/funeigen/utils/utils.R")

estimate_bandwidth_evalues <- function(curves, grid_bandwidth, grid_smooth, k0,
                                       grid_param, sigma = NULL, mu0 = NULL, 
                                       nvalues) 
  {
  #first use gkp covariance to estimate "pointwise" eigenfunctions
  cov_gkp <- covariance_ll(curves, grid_bandwidth, grid_smooth, k0,
                           grid_param, sigma, mu0)
  
  #obtain the normalised eigenvalues
  eigen_elements <- normalise_eigen(cov_gkp$Gamma, nvalues)
  
  #extract normalised eigenvalues
  evalues_norm <- eigen_elements$values
  efunctions_norm <- eigen_elements$vectors
  
  h_alpha <- sapply(grid_bandwidth, function(h) {
    cov_gkp$bw$params$constants$L * h^(2*cov_gkp$bw$params$constants$H)
  })
  
  q1_ts_t <- sapply(seq(nvalues), function(j) {
    efunctions_norm[, j]^2 * h_alpha * cov_gkp$bw$params$kernel_int
  }, simplify = "array")
  
  q1_ts_s_int <- sapply(seq(nvalues), function(j) {
    cov_gkp$bw$params$moments$EXt2 * abs(efunctions_norm[, j])
  }) |>
    apply(2, function(s) pracma::trapz(grid_smooth, s))

  q1_ts_t_int <- apply(q1_ts_t, c(2, 3),
                       function(t) pracma::trapz(grid_smooth, t))

  q1_ts_int <- sapply(seq_along(q1_ts_s_int), function(j) {
    q1_ts_t_int[, j] * q1_ts_s_int[j]
  })
  
  #q1_ts = q1_st here
  q1_term <- 4 * q1_ts_int
  
  q2_ts_t <- sapply(seq_along(grid_bandwidth), function(h) {
    sapply(seq(nvalues), function(j) {
      matrix(efunctions_norm[, j]^2, nrow = length(grid_smooth),
             ncol = length(grid_smooth)) * (1/cov_gkp$bw$params$Ngamma_ts)[,,h]
      }, simplify = "array")
    }, simplify = "array")
  
  #do everything in one shot - this gives everything before integrating
  #G x G x J x H
  q2_ts_int <- sapply(seq_along(grid_bandwidth), function(h) {
    sapply(seq(nvalues), function(j) {
      matrix(efunctions_norm[, j]^2, ncol = length(grid_smooth),
             nrow = length(grid_smooth)) *
        q2_ts_t[,,j,h] * matrix(cov_gkp$bw$params$moments$EXt2 * 
                                  efunctions_norm[, j]^2, 
                                ncol = length(grid_smooth),
                                nrow = length(grid_smooth),
                                byrow = TRUE)
    }, simplify = "array")
  }, simplify = "array")
  
  #to get s|t term, simply permute our arrays along the 1st and 2nd dimension
  #integrate out by t, and then s
  q2_term <- (q2_ts_int + aperm(q2_ts_int, c(2,1,3,4))) |>
    apply(c(2,3,4), function(x) pracma::trapz(grid_smooth, x)) |>
    apply(c(2,3), function(x) pracma::trapz(grid_smooth, x)) |>
    t()
  
  q3 <- array(cov_gkp$bw$params$moments$EXtXs2, 
                   dim = c(length(grid_smooth),
                           length(grid_smooth),
                           length(grid_bandwidth))) *
    ((1 / cov_gkp$bw$params$WN) - (1 / length(curves)))
  
  q3_term <- sapply(seq(nvalues), function(j) {
    sapply(seq_along(grid_bandwidth), function(h) {
      q3[,,h] * matrix(efunctions_norm[, j]^2, nrow = length(grid_smooth),
                       ncol = length(grid_smooth)) * 
        matrix(efunctions_norm[, j]^2, nrow = length(grid_smooth),
               ncol = length(grid_smooth), byrow = TRUE)
    }, simplify = "array")
  }, simplify = "array") |>
    apply(MARGIN = c(2,3,4), function(x) pracma::trapz(grid_smooth, x)) |>
    apply(MARGIN = c(2,3), function(x) pracma::trapz(grid_smooth, x))
    
  risk <- q1_term + q2_term + q3_term
  
  min_h_index <- apply(risk, MARGIN = 2, which.min)
  
  sapply(min_h_index, function(id) grid_bandwidth[id])
} 


smooth_curves_evalues <- function(curves, grid_bandwidth, grid_smooth, k0,
                                  grid_param, sigma = NULL, mu0 = NULL, 
                                  nvalues)
{
  
  bw_vector <- estimate_bandwidth_evalues(curves, grid_bandwidth, grid_smooth, k0,
                                          grid_param, sigma, mu0, nvalues)
  
  M_length <- curves |> purrr::map_dbl(~(length(.x$t)))
  
  #dim Mi x G x J for each curve i
  Wm_num <- purrr::map(curves, ~outer(.x$t, grid_smooth, FUN = "-")) |>
    purrr::map(~outer(.x, bw_vector, FUN = "/")) |>
    purrr::map(~epa_kernel(.x))
  
  Wm_denom <- purrr::map(Wm_num, ~colSums(.x, na.rm = TRUE))
  
  Wm_num_t <- purrr::map(Wm_num, ~aperm(.x, c(2, 1, 3)))
  
  Wm <- purrr::map2(Wm_num_t, Wm_denom, function(x, y) {
    sapply(seq(nvalues), function(j) {
      x[,,j] / y[, j]
    }, simplify = "array")
  })
  
  Wm <- lapply(Wm, function(w) replace(w, is.nan(w), 0))
  
  Xt <- purrr::map2(Wm, curves, function(w, y) {
    sapply(seq(nvalues), function(j) {
      w[,,j] %*% y$x
    })
  })
  list(smoothed_curves = Xt, bw = bw_vector)
}

evalues_adaptive <- function(curves, grid_bandwidth, grid_smooth, k0,
                             grid_param, sigma = NULL, mu0 = NULL, 
                             nvalues = 10) {
  
  smooth_curves <- smooth_curves_evalues(curves, grid_bandwidth,
                                         grid_smooth, k0,
                                         grid_param, sigma, mu0, nvalues)
  
  mu_eigen <- mean_plugin_evalues(curves, smooth_curves$smoothed_curves,
                                  smooth_curves$bw, grid_smooth, k0, nvalues)
  
  
  centered_curves <- lapply(smooth_curves$smoothed_curves, function(i) {
    i - mu_eigen$mu
  })
  
  weighted_curves <- purrr::map2(centered_curves, mu_eigen$wt,
                                 ~.x * .y)
  
  
  weighted_cov <- lapply(weighted_curves, function(i) {
    sapply(seq(nvalues), function(j) {
      i[, j] %*% t(i[, j])
    }, simplify = "array")
  }) |> (\(x) Reduce('+', x))()
  
  WN <- lapply(mu_eigen$wt, function(i) {
    sapply(seq(nvalues), function(j) {
      i[, j] %*% t(i[, j])
    }, simplify = "array")
  }) |>
    (\(x) Reduce('+', x))()
  
  emp_cov <- (1/WN) * weighted_cov
  
  
  eelements <- lapply(seq(nvalues), function(j) {
    normalise_eigen(emp_cov[,,j], nvalues)
  })
  
  evalues <- sapply(seq_along(eelements), function(j) {
    eelements[[j]]$values[j]
  })
  
  efunctions <- sapply(seq_along(eelements), function(j) {
    eelements[[j]]$vectors[,j]
  })
  
  list(eelements = list(values = evalues,
                        functions = efunctions),
      cov = emp_cov,
      bw = smooth_curves$bw,
      smooth_curves = smooth_curves$smoothed_curves)
}



