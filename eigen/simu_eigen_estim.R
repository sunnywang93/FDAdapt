#import dependent scripts
source("/Users/swang/Dropbox/funeigen/generation/generation_klutchnikoff.R")
source("/Users/swang/Dropbox/funeigen/regularity/H_and_L_functions.R")
source("/Users/swang/Dropbox/funeigen/mean/mean_optimised.R")
source("/Users/swang/Dropbox/funeigen/covariance/cov_optimised.R")
source("/Users/swang/Dropbox/funeigen/utils/utils.R")
source("/Users/swang/Dropbox/funeigen/eigen/eigenvalues.R")

library(dplyr)
library(purrr)
library(ggplot2)
library(reshape2)
library(foreach)
library(fs)
library(doParallel)
library(functional)
#parameters
N <- 50
M <- 30
H <- c(0.5, 0.8)
sigma <- 0.1
L <- 1
mu0 <- 1
grid_size_true <- 201
grid_true <- seq(0, 1, length.out = grid_size_true)
grid_mu <- seq(0, 1, length.out = 1001)
grid_smooth <- seq(0, 1, length.out = 101)
grid_bandwidth <- lseq(0.005, 0.15, length.out = 151)
grid_param = seq(0.1, 0.9, length.out = 20)
k0 <- 1
n_simu <- 1000
k_length <- 100
alpha <- 2
change_point <- 0.5
slope <- 50
tau <- 1
scale <- 1
shift <- 0
nvalues <- 10
#determine distribution of Ti - don't forget to specify additional parameters
points_dist <- Curry(runif, min = 0, max = 1)
set.seed(123)
mu_t <- generate_mean_curve(k_length = k_length, grid_t = grid_mu, alpha = alpha,
                            shift = shift, scale_mu = scale)
idx <- map_dbl(grid_smooth, ~which.min(abs(.x - mu_t$t)))
mu_Ti <- mu_t$mu[idx]
#no of sample paths you want to plot later 
NN <- 10 
#plot mean curve
plot(grid_smooth, mu_Ti, type = "l", xlab = "Ti", 
     main = "True Mean Curve (on smoothed grid)", col = "blue")



#determine the transition path of H - choose between "logistic" or "piecewise"
hurst_type <- "logistic"
if(hurst_type == "logistic") {
  hurst <- purrr::partial(hurst_logistic, h_left = !!H[1], h_right = !!H[2], 
                          change_point_position = !!change_point,
                          slope = !!slope)
} else {
  hurst <- purrr::partial(hurst_piecewise, h_left = !!H[1], h_right = !!H[2],
                          change_points_position = !!change_point)
}


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
    each_filename <- paste0('RESULT_N', N, '_M', M, '_', as.character(i),
                            '_', 'H(', H[1], ',', H[2],')', '.rda')
    each_filepath <- file.path(intermediate_directory, each_filename)
    
    if (file.exists(each_filepath)) {
      next
    }
    
    points_list <- generates_points(N, M, distribution = points_dist)
    
    pfbm <- generate_curves(points_list, hurst, 
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
    
    
    sample_seq <- runif(NN, min = 1, max = N) |> round() |> sort()
    
    result <- list(
      'pfbm_sampled' = pfbm$curves[sample_seq],
      't_min' = pfbm$t_min,
      't_max' = pfbm$t_max,
      'x_min' = pfbm$x_min,
      'x_max' = pfbm$x_max,
      'sample_seq' = sample_seq,
      'pfbm_mu' = pfbm_mu,
      'evalues_pw' = evalues_pw,
      'evalues_pointwise' = evalues_pointwise,
      'evalues_zhang' = evalues_zhang
      #'cov_gkp' = cov_gkp,
      #'cov_zhang' = cov_zhang,
      #'cov_emp' = cov_emp
    )
    save(result, file = each_filepath)
  }

tictoc::toc()
stopCluster(cl)

fls <- list.files(intermediate_directory, pattern = 'rda')
result_list <- lapply(fls,
                      function(x) get(eval(load(paste0(intermediate_directory, '/', x
                      )))))
save(result_list, file = paste0('evalues_', N, '_M', M, '_H(', H[1], ',', H[2],')', '.rda'))
path <- paste0(getwd(), '/intermediate')
file_delete(path)

#monte carlo simulation to obtain target eigenvalues
N_true <- 50000
M_true <- length(grid_smooth)

points_list_true <- generates_points(N_true, M_true, distribution = points_dist)
pfbm_true <- generate_curves(points_list_true, hurst, 
                             add_regular_grid_of_size = M_true,
                             sigma = 0,
                             tau = tau,
                             L = L)
pfbm_reshaped_true <- map(pfbm_true$curves, 
                          ~list(t = .x$ideal$t, x = .x$ideal$x))

#target covariance function
cov_true <- map(pfbm_reshaped_true, ~(.x$x %*% t(.x$x))) |>
  (\(x) Reduce('+', x) / (length(x) - 1))()

eigen_elements_true <- eigen(cov_true, symmetric = TRUE)

evalues_true_norm <- eigen_elements_true$values / M_true #(length_grid / M_true)
evectors_true_norm <- sapply(seq(M_true), function(k) {
  eigen_elements_true$vectors[, k] * sqrt(M_true) #sqrt(M_true / length_grid)
})


#analysis of results ===========================================
errors <- sapply(result_list, function(res) {
  res$evalues_pw - evalues_true_norm[seq(nvalues)]
})

errors_pointwise <- sapply(result_list, function(res) {
  res$evalues_pointwise - evalues_true_norm[seq(nvalues)]
})

errors_zw <- sapply(result_list, function(res) {
  res$evalues_zhang - evalues_true_norm[seq(nvalues)]
})

ratio_zw <- log(abs(errors) / abs(errors_zw))
ratio_pointwise <- log(abs(errors) / abs(errors_pointwise))

plot_list <- list()

for(j in 1:nvalues) {
  plot_list[[j]] <- tibble(ratio_zw = ratio_zw[j,], 
                           ratio_pointwise = ratio_pointwise[j,]) |>
    reshape2::melt() |>
    ggplot(aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, col = "red") +
    scale_y_continuous(breaks = seq(-5, 3, by = 0.5)) +
    xlab(paste0("eigenvalue", j)) +
    ylab("Log absolute ratio")
} 
    

plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,
                        guides = "collect") +
    plot_annotation(title = paste0("R = ", n_simu, ", N = ", N, ", M = ", M,
                                   ", H = ", H[1], ", ", H[2], ", scale = ", 
                                   scale, ", shift = ", shift, ", L2 = ", L^2))
}

plot_a_list(plot_list, 2, 5)





