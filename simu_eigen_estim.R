#import dependent scripts
#source("/Users/swang/Dropbox/funeigen/generation/generation_klutchnikoff.R")
source("/Users/swang/Dropbox/FPCA/simuls/nice_dgp/dgp.R")
# source("/Users/swang/Dropbox/funeigen/R/H_and_L_functions.R")
# source("/Users/swang/Dropbox/funeigen/R/mean_optimised.R")
# source("/Users/swang/Dropbox/funeigen/R/cov_optimised.R")
# source("/Users/swang/Dropbox/funeigen/R/utils.R")
# source("/Users/swang/Dropbox/funeigen/R/eigenvalues.R")

library(dplyr)
library(purrr)
library(ggplot2)
library(reshape2)
library(foreach)
library(fs)
library(doParallel)
library(functional)
#parameters
N <- 200
M <- 25
#H <- c(0.5, 0.8)
sigma <- 0.1
L <- 2
mu0 <- 1
grid_size_true <- 101
grid_true <- seq(0, 1, length.out = grid_size_true)
grid_mu <- seq(0, 1, length.out = 1001)
grid_smooth <- seq(0, 1, length.out = 101)
grid_bandwidth <- lseq(0.01, 0.15, length.out = 151)
grid_param = seq(0.1, 0.9, length.out = 20)
k0 <- 1
n_simu <- 45
# k_length <- 100
# alpha <- 2
#change_point <- 0.5
#slope <- 50
tau <- 1 #random starting points of parameters
# scale <- 1
# shift <- 0
nvalues <- 10
points_dist <- Curry(runif, min = 0, max = 1) #distribution of Ti

#learn true quantities from dataset powerconsumption ===========================
library(simulater)
library(fda)
df <- powerconsumption

df_list <- list()
for (idx in 1:nrow(df)) {
  df_list[[idx]] <- list(
    't' = seq(0, 1, length.out = 1440),
    'x' = unname(unlist(as.vector(df[idx,])))
  )
}

n_basis <- 9  #other tested setups 7
t0_list <- seq(.1, .9, l = 40)
grid <- seq(.1, .9, l = 100)

df_smooth <- funestim::presmoothing(df_list, t0_list = t0_list)
H0 <- funestim::estimate_H0(df_smooth)
L0 <- funestim::estimate_L0(df_smooth, H0, ncol(df))
variance <- funestim::estimate_var(df_smooth)

basis <- create.fourier.basis(rangeval = c(min(t0_list), max(t0_list)),
                              nbasis = n_basis)
H0_smooth <- smooth.basis(argvals = t0_list, y = H0, fdParobj = basis)
L0_smooth <- smooth.basis(argvals = t0_list, y = L0, fdParobj = basis)
H0_eval <- eval.fd(evalarg = grid, fdobj = as.fd(H0_smooth))[, 1]
L0_eval <- eval.fd(evalarg = grid, fdobj = as.fd(L0_smooth))[, 1]

hurst_fun <- approxfun(
  x = grid, y = H0_eval,
  yleft = H0_eval[1], yright = H0_eval[length(H0_eval)]
)
constant_fun <- approxfun(
  x = grid, y = L0_eval,
  yleft = L0_eval[1], yright = L0_eval[length(L0_eval)]
)

grid_large <- seq(0, 1, l = 10001)
A_hat_prime <- constant_fun(grid_large)**(1 / hurst_fun(grid_large))
A_hat <- pracma::cumtrapz(x = grid_large, y = A_hat_prime)
A_fun <- approxfun(
  x = grid_large, y = A_hat,
  yleft = A_hat[1], yright = A_hat[length(A_hat)]
)

s <- exp(-3)  # Regularity for the estimation of the mean
##other tested setups 3
hurst_fun <- approxfun(
  x = grid, y = H0_eval,
  yleft = H0_eval[1], yright = H0_eval[length(H0_eval)]
)

#true mean
mu_model <- learn_mean(df, k = 50)
mu <- predict_mean(grid_mu, mu_model, lambda = s, k = 50) - 240
mu_t <- tibble(t = grid_mu, mu)

#tested discretisation error on grid of 1001 points - very small
grid_cov <- seq(0, 1, length.out = 101)
pp_disto <- A_fun(grid_cov)
cov_true <- tau**2 + covariance_mfbm(grid_cov, hurst_fun, pp_disto)

#obtain the "numerically true" eigenvalues from the true covariance
evalues_true <- eigen(cov_true, symmetric = TRUE,
                      only.values = TRUE)$values[1:nvalues] / length(grid_cov)
efunctions_true <- eigen(cov_true, symmetric = TRUE)$vectors[, 1:nvalues] *
  sqrt(length(grid_cov))


#================================================================================


#### start simulations #########
intermediate_directory2 <- './intermediate'
if (!dir.exists(intermediate_directory2)){
  dir.create(intermediate_directory2)
}

cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)
tictoc::tic()

foreach::foreach(i = 1:n_simu,
                 .packages = c("dplyr", "MASS", "rpart", "matrixStats",
                               "purrr", "functional")) %dopar%
  {
    each_filename <- paste0('RESULT_N', N, '_M', M, '_', as.character(i),
                            '.rda')
    each_filepath <- file.path(intermediate_directory2, each_filename)

    if (file.exists(each_filepath)) {
      next
    }

    points_list <- generates_points(N, M, distribution = points_dist)

    pfbm <- generate_curves(points_list, hurst_fun,
                            distortion_model = A_fun,
                            add_regular_grid_of_size = grid_size_true,
                            sigma = sigma,
                            tau = tau,
                            L = L)


    pfbm_curves <- lapply(pfbm$curves,
                          function(x) list(t = x$observed$t ,x = x$observed$x))

    pfbm_mu <- add_mean_curve(data = pfbm_curves, mu_t = mu_t)

    evalues_pw <- evalues_adaptive(pfbm_mu, grid_bandwidth,
                                   grid_smooth, k0, grid_param,
                                   sigma = sigma, mu0 = mu0,
                                   nvalues = nvalues)

    cov_gkp <- covariance_ll(pfbm_mu, grid_bandwidth,
                             grid_smooth, k0, grid_param,
                             sigma = sigma, mu0 = mu0)

    evalues_pointwise <- normalise_eigen(cov_gkp$Gamma,
                                         nvalues = nvalues)$values

    cov_zhang <- covariance_lll(pfbm_mu, grid_smooth)

    evalues_zhang <- normalise_eigen(cov_zhang, nvalues = nvalues)$values

    evalues_intp <- eigen_interpolate(pfbm_mu, grid_smooth,
                                    nvalues = nvalues)$evalues

    evalues_kneip <- eigen_kneip(pfbm_mu, grid_smooth,
                               nvalues = nvalues)$evalues


    result <- list(
      'pfbm_mu' = pfbm_mu,
      'evalues_pw' = evalues_pw,
      'evalues_pointwise' = evalues_pointwise,
      'evalues_inter' = evalues_intp,
      'evalues_kneip' = evalues_kneip
    )
    save(result, file = each_filepath)
  }

tictoc::toc()
stopCluster(cl)

fls <- list.files(intermediate_directory2, pattern = 'rda')
result_list <- lapply(fls,
                      function(x) get(eval(load(paste0(intermediate_directory2, '/', x
                      )))))
save(result_list, file = paste0('evalues_', N, '_M', M, '.rda'))
path <- paste0(getwd(), '/intermediate')
file_delete(path)


#analysis of results ===========================================
errors <- sapply(result_list, function(res) {
  res$evalues_pw$eelements$values - evalues_avg
})

errors_emp <- sapply(result_list, function(res) {
  res$evalues_pw$eelements$values - res$eelements_emp$values
})

errors_pointwise <- sapply(result_list, function(res) {
  res$evalues_pointwise - evalues_avg
})

errors_pointwise_emp <- sapply(result_list, function(res) {
  res$evalues_pointwise - res$eelements_emp$values
})

errors_zw <- sapply(result_list, function(res) {
  res$evalues_zhang - evalues_avg
})

errors_zw_emp <- sapply(result_list, function(res) {
  res$evalues_zhang - res$eelements_emp$values
})

errors_emp_true <- sapply(result_list, function(j) {
  j$eelements_emp$values - evalues_avg
})

#ratio_zw <- log(abs(errors) / abs(errors_zw))
#ratio_pointwise <- log(abs(errors) / abs(errors_pointwise))
ratio_zw <- log(abs(errors) / abs(errors_zw))
ratio_pointwise <- log(abs(errors) / abs(errors_pointwise))

ratio_zw_emp <- log(abs(errors_emp) / abs(errors_zw_emp))
ratio_pointwise_emp <- log(abs(errors_emp) / abs(errors_pointwise_emp))

ratio_emp_true <- log(abs(errors) / abs(errors_emp_true))

plot_list <- list()

for(j in 1:nvalues) {
  plot_list[[j]] <- tibble(ratio_zw = ratio_zw[j,],
                           ratio_pointwise = ratio_pointwise[j,]) |>#,
                           #ratio_zw_emp = ratio_zw_emp[j,],
                           #ratio_pointwise_emp = ratio_pointwise_emp[j,],
                           #ratio_emp_true = ratio_emp_true[j,]) |>
    reshape2::melt() |>
    ggplot(aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, col = "red") +
    scale_y_continuous(breaks = seq(-5, 3, by = 0.5)) +
    xlab(paste0("eigenvalue", j)) +
    ylab("Log Ratio Absolute Error") +
    theme(legend.position = "none")
}


plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {

  patchwork::wrap_plots(master_list_with_plots,
                        nrow = no_of_rows, ncol = no_of_cols,
                        guides = "collect") +
    patchwork::plot_annotation(title = paste0("R = ", n_simu, ", N = ", N, ", M = ", M,
                                   ", H = (", H[1], ",", H[2], "),", " L = ", L,
                                   ", tau = ", tau, ", sigma = ", sigma))
}

plot_a_list(plot_list, 2, 5)

