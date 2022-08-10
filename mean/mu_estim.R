source("/Users/swang/Dropbox/research_code/generation_klutchnikoff.R")
source("/Users/swang/Dropbox/research_code/H_and_L_functions.R")
source("/Users/swang/Dropbox/research_code/mean_optimised.R")
source("/Users/swang/Dropbox/research_code/utils.R")
library(dplyr)
library(purrr)
library(ggplot2)
library(reshape2)
library(foreach)
library(fs)
library(doParallel)
library(functional)
#parameters
N <- 100
M <- 100
H <- c(0.5, 0.8)
sigma <- 0.1
L <- c(4, 1)
mu0 <- 1
grid_size_true <- 201
grid_mu <- seq(0, 1, length.out = 1001)
grid_smooth <- seq(0, 1, length.out = 101)
grid_bandwidth <- lseq(0.005, 0.1, length.out = 151)
k0 <- 2
n_simu <- 50
k_length <- 60
alpha <- 2
change_point <- 0.5
slope <- 50
tau <- 1
#determine distribution of Ti - don't forget to specify additional parameters
#if required (e.g shape in beta distribution): pass these as args to Curry
points_dist <- Curry(runif, min = 0, max = 1)
set.seed(123)
mu_t <- generate_mean_curve(k_length = k_length, grid_t = grid_mu, alpha = alpha)
mu_Ti <- map_dbl(grid_smooth, ~which.min(abs(.x - mu_t$t)))
mu_Ti <- mu_t$mu[mu_Ti]
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

#start simulation
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
                            tau = tau) 
    
    
    pfbm_curves <- lapply(pfbm$curves, 
                   function(x) list(t = x$observed$t ,x = x$observed$x))

    pfbm_mu <- add_mean_curve(data = pfbm_curves, mu_t = mu_t)
    holder_const <- estimate_holder_const(pfbm_mu, sigma = sigma, mu0 = mu0)
    
    c <- Sys.time()
    mu_gkp <- mean_ll(pfbm_mu, holder_const, grid_bandwidth = grid_bandwidth,
                      grid_smooth = grid_smooth, k0 = k0)
    time_gkp <- difftime(Sys.time(), c, units = "secs")
    
    c <- Sys.time()
    mu_cai <- mean_ss(pfbm_mu, grid_smooth)
    time_cai <- difftime(Sys.time(), c, units = "secs")
    c <- Sys.time()
    mu_zhang <- mean_lll(pfbm_mu, grid_smooth)
    time_zhang <- difftime(Sys.time(), c, units = "secs")

    
    result <- list(
      'pfbm_all' = pfbm,
      'pfbm_mu' = pfbm_mu,
      'holder_const' = holder_const,
      'mu_gkp' = mu_gkp,
      'mu_cai' = mu_cai,
      'mu_zhang' = mu_zhang,
      'time_gkp' = time_gkp,
      'time_cai' = time_cai,
      'time_zhang' = time_zhang
    )
    save(result, file = each_filepath)
  }

tictoc::toc()
stopCluster(cl)

fls <- list.files(intermediate_directory, pattern = 'rda')
result_list <- lapply(fls,
                      function(x) get(eval(load(paste0(intermediate_directory, '/', x
                      )))))
save(result_list, file = paste0('mu_', N, '_M', M, '_H(', H[1], ',', H[2],')', '.rda'))
path <- paste0(getwd(), '/intermediate')
file_delete(path)



#plot one random sample (of simulation) of the mean function
sample_id <- runif(1, min = 1, max = n_simu) |> round()
tibble(t = grid_smooth, gkp = result_list[[sample_id]]$mu_gkp$mu_hat,
       cai = result_list[[sample_id]]$mu_cai,
       zhang = result_list[[sample_id]]$mu_zhang,
       true = mu_Ti) |>
  melt(id = "t") |> ggplot(aes(x = t, y = value, col = variable)) +
  geom_line()

#plot some sample paths of the generated process
NN <- 10 #change this if you want to plot more curves

plot(result_list[[sample_id]]$pfbm_all,
  which_curves = 1:NN,
  main =  paste(NN, "Sample Paths ( simu =", sample_id, ")"),
  xlab = "time",
  ylab = "value")



#edit below according to N and M to get your desired result 
mean_results <- get_means(result_list, n_simu = n_simu)
tibble(ratio_CY = mean_results[, 1], ratio_ZW = mean_results[, 2]) %>%
  reshape2::melt() %>%
  ggplot(aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  labs(title = paste0("R = ", n_simu, ", N = ", N, ", M = ", M,
                      ", alpha = ", alpha))






