#import dependent scripts
source("/Users/swang/Dropbox/funeigen/generation/generation_klutchnikoff.R")
source("/Users/swang/Dropbox/funeigen/regularity/H_and_L_functions.R")
source("/Users/swang/Dropbox/funeigen/mean/mean_optimised.R")
source("/Users/swang/Dropbox/funeigen/covariance/cov_optimised.R")
source("/Users/swang/Dropbox/funeigen/utils/utils.R")

#load dependent packages
library(dplyr)
library(purrr)
library(ggplot2)
library(reshape2)
library(foreach)
library(fs)
library(doParallel)
library(functional)

#monte-carlo simulation to obtain target eigenvalues ===========================

#parameters
grid_smooth <- seq(0, 1, length.out = 101)
#length_grid <- max(grid_smooth) - min(grid_smooth)
k_length <- 100
grid_mu <- seq(0, 1, length.out = length(grid_smooth))
alpha <- 2
shift <- 0
scale <- 1
N_true <- 2000
M_true <- length(grid_smooth)
points_dist <- Curry(runif, min = 0, max = 1)
H <- c(0.5, 0.8)
tau <- 1
L <- 2
sigma <- 0.1
change_point <- 0.5
slope <- 50
N <- 50
M <- 30
grid_bandwidth <- lseq(0.005, 0.15, length.out = 151)
k0 <- 1
grid_param <- seq(0.1, 0.9, length.out = 20)
mu0 <- 1
NN <- 10 #no of sample paths to plot in each simulation
n_simu <- 50


#generate true mean curve
set.seed(123)
mu_t <- generate_mean_curve(k_length = k_length, grid_t = grid_mu, alpha = alpha,
                            shift = shift, scale_mu = scale)

#no of sample paths to plot later 
NN <- 10 
#plot mean curve
plot(mu_t$t, mu_t$mu, type = "l", xlab = "t", ylab = "x",
     main = "True Mean Curve", col = "blue")

#logistic change for regularity
hurst <- purrr::partial(hurst_logistic, h_left = !!H[1], h_right = !!H[2], 
                        change_point_position = !!change_point,
                        slope = !!slope)

#generate target curves - noiseless on common design, takes some time
#due to large sample size
points_list_true <- generates_points(N_true, M_true, distribution = points_dist)
pfbm_true <- generate_curves(points_list_true, hurst, 
                        add_regular_grid_of_size = M_true,
                        sigma = 0,
                        tau = tau,
                        L = L)
pfbm_reshaped_true <- map(pfbm_true$curves, ~list(t = .x$ideal$t, x = .x$ideal$x))

#target covariance function
cov_true <- map(pfbm_reshaped_true, ~(.x$x %*% t(.x$x))) |>
  (\(x) Reduce('+', x) / (length(x) - 1))()

eigen_elements_true <- eigen(cov_true, symmetric = TRUE)

evalues_true_norm <- eigen_elements_true$values / M_true #(length_grid / M_true)
evectors_true_norm <- sapply(seq(M_true), function(k) {
  eigen_elements_true$vectors[, k] * sqrt(M_true) #sqrt(M_true / length_grid)
})
  
#run simulations for estimation of curves ======================================

#to ensure we have the same seed, we can draw a vector of seeds first and then
#use those in the for loop for each run of the simulation below

intermediate_directory <- './intermediate'
if (!dir.exists(intermediate_directory)){
  dir.create(intermediate_directory)
}

cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
tictoc::tic()

foreach::foreach(i = 1:n_simu, 
                 .packages = c("dplyr", "MASS", "rpart", "matrixStats",
                               "purrr", "functional")) %dopar%
  {
    each_filename <- paste0('EIGENRESULT_N', N, '_M', M, '_', as.character(i),
                            '_', 'H(', H[1], ',', H[2],')', '.rda')
    each_filepath <- file.path(intermediate_directory, each_filename)
    
    if (file.exists(each_filepath)) {
      next
    }
    
    points_list <- generates_points(N, M, distribution = points_dist)
    pfbm <- generate_curves(points_list, hurst, 
                            add_regular_grid_of_size = M_true,
                            sigma = sigma,
                            tau = tau,
                            L = L)
    pfbm_reshaped <- lapply(pfbm$curves, 
                          function(x) list(t = x$observed$t ,x = x$observed$x))
    pfbm_reshaped_mean <- add_mean_curve(data = pfbm_reshaped, mu_t = mu_t)
    
    #estimate eigenelements with different covariance estimators
    cov_gkp <- covariance_ll(pfbm_reshaped_mean, grid_bandwidth,
                             grid_smooth, k0,
                             grid_param, sigma, mu0)
    eigen_gkp <- eigen(cov_gkp$Gamma, symmetric = TRUE)
    
    cov_zhang <- covariance_lll(pfbm_reshaped_mean, grid_smooth)
    eigen_zhang <- eigen(cov_zhang, symmetric = TRUE)
    
    sample_seq <- runif(NN, min = 1, max = N) |> round() |> sort()
    
    result <- list(
      'pfbm_sampled' = pfbm$curves[sample_seq],
      'cov_gkp' = cov_gkp,
      'cov_zhang' = cov_zhang,
      'eigen_gkp' = eigen_gkp,
      'eigen_zhang' = eigen_zhang
    )
    save(result, file = each_filepath)
  }

tictoc::toc()
stopCluster(cl)

fls <- list.files(intermediate_directory, pattern = 'rda')
result_list <- lapply(fls,
                      function(x) get(eval(load(paste0(intermediate_directory, '/', x
                      )))))
save(result_list, file = paste0('EIGEN_', N, '_M', M, 
                                '_H(', H[1], ',', H[2],')', '.rda'))
path <- paste0(getwd(), '/intermediate')
file_delete(path)


#analysis of simulation results ================================================

#extract results
evalues_gkp <- sapply(result_list, function(x) x$eigen_gkp$values)
evectors_gkp <- sapply(result_list, function(x) x$eigen_gkp$vectors,
                       simplify = "array")

evalues_zw <- sapply(result_list, function(x) x$eigen_zhang$values)
evectors_zw <- sapply(result_list, function(x) x$eigen_zhang$vectors,
                      simplify = "array")

#normalise eigenvalues
evalues_gkp_norm <- sapply(seq(n_simu), 
                           function(i) evalues_gkp[, i] / length(grid_smooth))

evalues_zw_norm <- sapply(seq(n_simu),
                          function(i) evalues_zw[, i] / length(grid_smooth))

#compare one eigenvalue
kj <- 1 #j-th eigenvalue to analyse
first_evalue_gkp_norm <- sapply(seq(n_simu),
                                function(i) evalues_gkp_norm[, i][kj])

first_evalues_zw_norm <- sapply(seq(n_simu),
                                function(i) evalues_zw_norm[, i][kj])

diff_firstevalue_gkp <- abs(first_evalue_gkp_norm - evalues_true_norm[kj])
diff_firstevalue_zw <- abs(first_evalues_zw_norm - evalues_true_norm[kj])

ratio <- log(diff_firstevalue_gkp /diff_firstevalue_zw)

tibble(ratio = ratio) |>
  melt() |>
  ggplot(aes(x = variable, y = value, fill = variable)) +
  geom_boxplot()


#normalise eigenvectors
evectors_gkp_norm <- sapply(seq(n_simu), 
                            function(i) evectors_gkp[,,i] * sqrt(M_true),
                            simplify = "array")


evectors_zw_norm <- sapply(seq(n_simu),
                           function(i) evectors_zw[,,i] * sqrt(M_true),
                           simplify = "array")

psikj <- 1 #j-th eigenvector to analyse
j_evector_gkp_norm <- sapply(seq(n_simu), 
                                 function(i) evectors_gkp_norm[,psikj,i])

sgn_gkp <- sapply(seq(n_simu), function(i) {
  c(sign(t(evectors_true_norm[, psikj]) %*% j_evector_gkp_norm[,i]))
  })

j_evector_zw_norm <- sapply(seq(n_simu), 
                                function(i) evectors_zw_norm[,psikj,i])


sgn_zw <- sapply(seq(n_simu), function(i) {
  c(sign(t(evectors_true_norm[, psikj]) %*% j_evector_zw_norm[,i]))
})

norm_gkp <- sapply(seq(n_simu), 
                   function(i) norm((sgn_gkp[i] * evectors_true_norm[, 1] - 
                          j_evector_gkp_norm[, i]) / sqrt(M_true), 
                          type = "2"))

norm_zw <- sapply(seq(n_simu), 
                  function(i) norm((sgn_zw[i] * evectors_true_norm[, 1] - 
                                      j_evector_zw_norm[, i]) / sqrt(M_true), 
                                   type = "2"))

ratio_j_evector <- log(norm_gkp / norm_zw)

tibble(ratio_j_evector) |>
  melt() |>
  ggplot(aes(x = variable, y = value, fill = variable)) +
  geom_boxplot()

