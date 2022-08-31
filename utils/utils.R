MSE_1D <- function(mu, mu_estim) {
  sum((mu - mu_estim)^2) / length(mu)
}

ISE <- function(t, x, y, trimmed = FALSE){
  if (trimmed) {
    cond <- (t > 0.05) & (t < 0.95)
    t <- t[cond]
    x <- x[cond]
    y <- y[cond]
  }
  if (length(t) == length(x) & length(x) == length(y)) {
    return(pracma::trapz(t, (x - y)**2))
  } else {
    return(NA)
  }
}

ISE_2D <- function(t, x, y, trimmed = FALSE) {
  x <- as.vector(x)
  y <- as.vector(y)
  if (trimmed) {
    cond <- (t > 0.05) & (t < 0.95)
    x <- x[cond]
    y <- y[cond]
  }
  if (length(x) == length(y)) {
    return(mean((x-y)**2))
  }
  else{
    return(NA)
  }
}

list2cai <- function(data){
  seq_along(data) |> 
    lapply(function(idx) {
      data.frame(obs = idx, time = data[[idx]]$t, x = data[[idx]]$x)
    }) |> 
    (\(x) do.call("rbind", x))()
}

mean_ss <- function(curves, grid){
  curves_ <- list2cai(curves)
  mod <- stats::smooth.spline(curves_$time, curves_$x)
  stats::predict(mod, grid)$y
}

mean_lll <- function(curves, grid) {
  curves_ <- list2cai(curves)
  L3 <- fdapace::MakeFPCAInputs(
    IDs = curves_$obs, tVec = curves_$time, yVec = curves_$x, deduplicate = TRUE,
    sort = TRUE)
  fdapace::GetMeanCurve(
    L3$Ly, L3$Lt,
    list(kernel = 'epan', nRegGrid = length(grid), methodBwMu = 'GCV')
  )$mu
}

generate_mean_curve <- function(k_length, grid_t = seq(0, 1, length.out = 101),
                                alpha, shift = 0, scale_mu = 1) {
  Z <- runif(k_length, -sqrt(3), sqrt(3))
  #Z <- rnorm(k_length)
  xi <- c()
  # for(i in 1:k_length) {
  #   xi_k <- (-1)^(i+1)*i^(-alpha)
  #   xi <- c(xi, xi_k)
  # }
  for(i in 1:k_length) {
    xi_k <- 1/((i - 1/2)^alpha*pi^2)
    xi <- c(xi, xi_k)
  }
  k <- seq(1, k_length, by = 1)
  
  # mu_kt <- sapply(grid_t, function(t) {
  #   Z * sqrt(xi) * cos(k * pi * t)
  # })
  
  mu_kt <- sapply(grid_t, function(t) {
    sqrt(2) * Z * sqrt(xi) * sin((k - 0.5) * pi * t)
  })
  

  mu_t <- scale_mu * colSums(mu_kt) + rep(shift, length(grid))
  tibble(t = grid_t, mu = mu_t)
} 


#add mean curve to each curve
#input: mu_t on a dense grid
add_mean_curve <- function(data, mu_t) {
  
  m <- data |> map_dbl(~length(.x$t))
  
  idx <- lapply(data, function(i) {
    map_dbl(i$t, ~which.min(abs(.x - mu_t$t)))
  })
  
  mu_Ti <- lapply(idx, function(id) {
    mu_t$mu[id]
  })
  
  curves_with_mean <- map2(data, mu_Ti, ~(.x$x + .y))
  
  for(i in 1:length(data)) {
    data[[i]]$x <- curves_with_mean[[i]]
  }
  data
}


add_mean_to_true <- function(data, mu_t) {
  idx <- purrr::map_dbl(data[[1]]$ideal$t, ~which.min(abs(.x - mu_t$t)))
  mu_grid_true <- mu_t$mu[idx]
  curves_with_mu <- purrr::map(data, ~.x$ideal$x + mu_grid_true)
  init <- vector(mode = "list", length = length(data))
  for (i in 1:length(init)) {
    init[[i]]$t <- data[[1]]$ideal$t
    init[[i]]$x <- curves_with_mu[[i]]
  }
  init
}


get_means <- function(result_list, n_simu = n_simu) {
  mu_gkp <- sapply(result_list, function(x) x$mu_gkp$mu_hat)
  mu_cai <- sapply(result_list, function(x) x$mu_cai)
  mu_zhang <- sapply(result_list, function(x) x$mu_zhang)
  
  MSE_gkp <- sapply(seq(n_simu), function(id) {
    MSE_1D(mu_Ti, mu_gkp[, id])
  })
  
  MSE_cai <- sapply(seq(n_simu), function(id) {
    MSE_1D(mu_Ti, mu_cai[, id])
  })
  
  MSE_zhang <- sapply(seq(n_simu), function(id) {
    MSE_1D(mu_Ti, mu_zhang[, id])
  })
  
  #ratios of integrated squared error
  ratio_cai <- MSE_gkp / MSE_cai
  ratio_zhang <- MSE_gkp / MSE_zhang
  
  cbind(ratio_cai, ratio_zhang)
}


covariance_lll <- function(curves, grid, bandwidth = 0.1){
  if (!inherits(curves, 'list')) curves <- checkData(curves)
  curves_ <- list2cai(curves)
  L3 <- fdapace::MakeFPCAInputs(
    IDs = curves_$obs, tVec = curves_$time, yVec = curves_$x,
    deduplicate = TRUE)
  fdapace::GetCovSurface(
    L3$Ly, L3$Lt,
    list(kernel = 'epan', nRegGrid = length(grid),
         methodMuCovEst = 'smooth', userBwCov = bandwidth)
  )$cov
}


true_cov_pfbm <- function(t, s, H, change_point) {
  obs <- expand.grid(t = t, s = s)
  0.5 * (abs(obs$t)^(2*H) + abs(obs$s)^(2*H) - abs(obs$t - obs$s)^(2*H))
}

normalise_eigen <- function(covariance, nvalues = 10) {
  eelements <- eigen(covariance, symmetric = TRUE)
  evalues <- eelements$values[seq(nvalues)]
  efunctions <- eelements$vectors[, seq(nvalues)]
  evalues_norm <- evalues / nrow(covariance)
  efunctions_norm <- sapply(seq(nvalues), function(j) efunctions[, j] * 
                              sqrt(nrow(covariance)))
  list(values = evalues_norm,
       vectors = efunctions_norm)
}

eigen_error <- function(eigen_estim, eigen_true) {
  evalues_error <- (eigen_estim$values - eigen_true$values) / eigen_true$values
  sgn <- sapply(seq(ncol(eigen_estim$vectors)), function(j) {
    c(sign(t(eigen_estim$vectors[, j]) %*% eigen_true$vectors[, j])) 
  })
  evectors_error <- sapply(seq(ncol(eigen_estim$vectors)), function(j) {
    norm((sgn[j] * eigen_true$vectors[, j] - eigen_estim$vectors[, j]) / 
           sqrt(nrow(eigen_estim$vectors)), type = "2") 
  })
  list(evalues_error = evalues_error, 
       efunctions_error = evectors_error)
}








