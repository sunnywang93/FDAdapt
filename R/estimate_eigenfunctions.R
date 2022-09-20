
estimate_bandwidth_efunctions <- function(curves, grid_bandwidth, grid_smooth,
                                          k0, grid_param, sigma, mu0,
                                          nfunctions)
  {

  cov_gkp <- covariance_ll(curves, grid_bandwidth, grid_smooth, k0,
                           grid_param, sigma, mu0)

  eigen_elements <- normalise_eigen(cov_gkp$Gamma, nvalues)

  evalues_norm <- eigen_elements$values
  efunctions_norm <- eigen_elements$vectors

  h_alpha <- sapply(grid_bandwidth, function(h) {
    cov_gkp$bw$params$constants$L * h^(2*cov_gkp$bw$params$constants$H)
  })




}






















