#import dependent scripts
source("/Users/swang/Dropbox/funeigen/generation/generation_klutchnikoff.R")
source("/Users/swang/Dropbox/funeigen/regularity/H_and_L_functions.R")
source("/Users/swang/Dropbox/funeigen/mean/mean_optimised.R")
source("/Users/swang/Dropbox/funeigen/covariance/cov_optimised.R")
source("/Users/swang/Dropbox/funeigen/utils/utils.R")


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
grid_bandwidth <- lseq(0.005, 0.1, length.out = 151)
k0 <- 2
n_simu <- 50
k_length <- 100
alpha <- 2
change_point <- 0.5
slope <- 50
tau <- 1
scale <- 1
shift <- 0
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
                               "purrr", "FKSUM", "functional")) %dopar%
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
    
    holder_const <- estimate_holder_const(pfbm_mu, sigma = sigma, mu0 = mu0)
    
    #compute empirical objects
    pfbm_mu_orig <- add_mean_to_original(pfbm$curves, mu_t, grid_smooth)
    
    mu_emp <- sapply(pfbm_mu_orig, function(i) i$x) |> rowMeans() 
    cov_emp <- purrr::map(pfbm_mu_orig, ~(.x$x - mu_emp) %*% (t(.x$x - mu_emp))) |>
      (\(x) Reduce('+', x) / (length(x) - 1))()
    
    cov_gkp <- covariance_ll(pfbm_mu, holder_const, grid_bandwidth,
                             grid_smooth, k0)
    
    cov_zhang <- covariance_lll(pfbm_mu, grid_smooth)

    
    sample_seq <- runif(NN, min = 1, max = N) |> round() |> sort()
    
    result <- list(
      'pfbm_sampled' = pfbm$curves[sample_seq],
      't_min' = pfbm$t_min,
      't_max' = pfbm$t_max,
      'x_min' = pfbm$x_min,
      'x_max' = pfbm$x_max,
      'sample_seq' = sample_seq,
      'pfbm_mu' = pfbm_mu,
      'holder_const' = holder_const,
      'cov_gkp' = cov_gkp,
      'cov_zhang' = cov_zhang,
      'cov_emp' = cov_emp
    )
    save(result, file = each_filepath)
  }

tictoc::toc()
stopCluster(cl)

fls <- list.files(intermediate_directory, pattern = 'rda')
result_list <- lapply(fls,
                      function(x) get(eval(load(paste0(intermediate_directory, '/', x
                      )))))
save(result_list, file = paste0('cov_', N, '_M', M, '_H(', H[1], ',', H[2],')', '.rda'))
path <- paste0(getwd(), '/intermediate')
file_delete(path)



#plot one random sample (of simulation) of the mean function
sample_id <- runif(1, min = 1, max = n_simu) |> round()


#plot some sample paths of the generated process
plot_subset(result_list[[sample_id]]$pfbm_sampled,
            t_min = result_list[[sample_id]]$t_min,
            t_max = result_list[[sample_id]]$t_max,
            x_min = result_list[[sample_id]]$x_min,
            x_max = result_list[[sample_id]]$x_max,
            main =  paste(NN, "Sample Paths ( simu =", sample_id, ")"),
            xlab = "time",
            ylab = "value")



#compute integrated squared errors
cov_true <- covariance_mfbm(grid_smooth, hurst) + tau^2 * 1


ISE_gkp_true <- purrr::map_dbl(result_list,
                               ~ISE_2D(grid_smooth,
                                .x$cov_gkp$Gamma,
                                 cov_true))

ISE_gkp_emp <- purrr::map_dbl(result_list,
                              ~ISE_2D(grid_smooth,
                               .x$cov_gkp$Gamma,
                               .x$cov_emp))

ISE_zhang_true <- purrr::map_dbl(result_list, 
                            ~ISE_2D(grid_smooth, 
                             .x$cov_zhang, cov_true))

ISE_zhang_emp <- purrr::map_dbl(result_list,
                                ~ISE_2D(grid_smooth,
                                 .x$cov_zhang, .x$cov_emp))


ISE_ratio_gkp_true <- ISE_gkp_true / ISE_zhang_true
ISE_ratio_gkp_emp <- ISE_gkp_emp / ISE_zhang_emp


emp <- tibble(ISE_ratio_gkp_true = ISE_ratio_gkp_emp) |>
  reshape2::melt() |>
  ggplot(aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, col = "red") +
  scale_y_continuous(breaks = seq(0, 5, by = 0.5)) +
  labs(title = "Empirical")


true <- tibble(ISE_ratio_gkp_true = ISE_ratio_gkp_true) |>
  reshape2::melt() |>
  ggplot(aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, col = "red") +
  scale_y_continuous(breaks = seq(0, 5, by = 0.5)) +
  labs(title = "True")

require("patchwork")
emp + true + plot_layout(guides = "collect") + 
  plot_annotation(title = paste0("R = ", n_simu, ", N = ", N, ", M = ", M,
                                 ", H = ", H[1], ", ", H[2], ", scale = ", scale,
                                 ", shift = ", shift, ", L2 = ", L^2))


require("GA")
par(mfrow = c(2, 2))
persp3D(grid_smooth, grid_smooth, cov_true, 
        xlab = "t", ylab = "s", zlab = "values",
        main = paste0("True Covariance (N = ", N, ", M = ", M, ")"),
        phi = 30,
        expand = 0.5)


persp3D(grid_smooth, grid_smooth, result_list[[sample_id]]$cov_emp, 
        xlab = "t", ylab = "s", zlab = "values",
        main = paste0("Empirical Covariance (N = ", N, ", M = ", M, ")"),
        phi = 30,
        expand = 0.5)


persp3D(grid_smooth, grid_smooth, result_list[[sample_id]]$cov_zhang, 
        xlab = "t", ylab = "s", zlab = "values",
        main = paste0("Zhang & Wang Covariance", " id = ", sample_id),
        phi = 30,
        expand = 0.5)


persp3D(grid_smooth, grid_smooth, result_list[[sample_id]]$cov_gkp$Gamma, 
        xlab = "t", ylab = "s", zlab = "values",
        main = paste0("GKP covariance", " id = ", sample_id),
        phi = 30,
        expand = 0.5)







