source("/Users/swang/Dropbox/research_code/H_and_L_functions.R")
library(dplyr)
library(purrr)
library(ggplot2)
library(reshape2)

lseq <- function (from = 1, to = 100, length.out = 51) {
  exp(seq(log(from), log(to), length.out = length.out))
}


neighbours <- function(t, t0, h, k0) {
  indic_mat <- sapply(t0, function(x) ifelse(abs(t - x) <= h, 1, 0))
  as.numeric(colSums(indic_mat) >= k0)
}


epa_kernel <- function(y) {
  3/4 * (1 - y^2) * (abs(y) <= 1)
}

laplace_kernel <- function(x) {
  #1/4 * exp(-abs(x)) + 1/4 * abs(x) * exp(-abs(x))
  exp(-abs(x)) + abs(x) * exp(-abs(x)) + 0.5 * abs(x)^2 * exp(-abs(x)) +
    1/6 * abs(x)^3 * exp(-abs(x)) + 1/24 * abs(x)^4 * exp(-abs(x))
}

bertin_kernel <- function(x, beta){
  ((1 + beta)/ (2 * beta))  * (1 - abs(x)^beta) * (abs(x) <= 1)
}

#input: params simply needs to include H and L
estimate_variance_curves <- function(data, params, grid_smooth,
                                     sigma = NULL, mu0 = NULL) {
  if(is.null(sigma)) {
    sigma <- estimate_sigma(data)
  }
  
  if(is.null(mu0)) {
    mu0 <- estimate_mu0(data)
  }
  
  m <- data |> purrr::map_dbl(~length(.x$t)) |> mean()
  #associate H on its own grid to the smoothing grid
  break_points <- params |> dplyr::group_by(H, L, sigma) |> 
    dplyr::summarise(t_break = max(t), .groups = "keep") |> 
    dplyr::arrange(t_break)
  
  intervals <- cut(grid_smooth, breaks = c(0, break_points$t_break, 1),
                   include.lowest = TRUE)
  intervals_H <- intervals
  levels(intervals_H) <- c(break_points$H, break_points$H[length(break_points$H)])
  
  intervals_L <- intervals
  levels(intervals_L) <- c(break_points$L, break_points$L[length(break_points$L)])
  
  intervals_sigma <- intervals
  levels(intervals_sigma) <- c(break_points$sigma, 
                               break_points$sigma[length(break_points$sigma)])
  grid_tibble <- tibble::tibble(t = grid_smooth, 
                                H = as.numeric(as.character(intervals_H)),
                                L = as.numeric(as.character(intervals_L)),
                                sigma = as.numeric(as.character(intervals_sigma)))
  #now smooth each curve using Bertin's bandwidth
  c <- (sigma^(2*grid_tibble$H) * grid_tibble$L * 
          ((grid_tibble$H + 1)/ 2 * grid_tibble$H^2 * mu0)^grid_tibble$H)^
    (1 / (2*grid_tibble$H + 1))
  
  psi_m <- (1 / m)^(grid_tibble$H / (2 * grid_tibble$H + 1))
  delta <- min(log(m)^(-1.1), 0.2)
  b_naive <- pmax(pmin((c * psi_m / grid_tibble$L)^(1 / grid_tibble$H),
                       delta/4), log(m)/m)
  
  
  
  #define nadaraya-watson weights 
  #here we use Bertin's bandwidth and kernel
  Wmi <- lapply(data, function(i) {
    sapply(i$t, function(Tmi) {
      bertin_kernel((Tmi - grid_smooth) / b_naive, beta = grid_tibble$H)
    }) |> (\(x) x / rowSums(x, na.rm = TRUE))() |>
      (\(x) {x[is.nan(x)] <- 0; x})()
  })

  
  #obtain smoothed curves 
  #output: G x N matrix
  X_hat <- mapply(function(W, Y) {
    W %*% Y$x
  }, Wmi, data)
  
  #might be a problem with numerical accuracy when computing manually
  #X_bar <- rowMeans(X_hat, na.rm = TRUE)
  # diff_var <- apply(X_hat, 2, function(x) (x - X_bar)^2)
  # var_Xt <- rowMeans(diff_var, na.rm = TRUE)
  var_Xt <- apply(X_hat, 1, var, na.rm = TRUE)
  E_Xt2 <- rowMeans(X_hat^2, na.rm = TRUE)
  #we have G X G matrix for each curve 
  X_hat_prod <- lapply(seq_along(data), function(i) {
    X_hat[, i] %*% t(X_hat[, i])
  })
  #use modifiedSum to ensure NAs do not affect the computation
  modifiedSum <- function(x, y) {
    replace(x, is.na(x), 0) + replace(y, is.na(y), 0)
  }
  
  X_bar_prod <- Reduce(modifiedSum, X_hat_prod)/length(X_hat_prod)
  
  diff_var_prod <- lapply(X_hat_prod, function(i) {
    (i - X_bar_prod)^2
  })
  
  var_XtXs <- Reduce(modifiedSum, diff_var_prod)/length(diff_var_prod)
  list(varXt = var_Xt, varXtXs = var_XtXs, EXt2 = E_Xt2)
}


#input: params must include H, L, sigma, mu0
estimate_bandwidth_mean <- function(curves, params, grid_bandwidth, 
                                    grid_smooth,
                                    k0) {
  
  #extrapolate values from params grid to grid_smooth
  break_points <- params |> dplyr::group_by(H, L, sigma) |> 
    dplyr::summarise(t_break = max(t), .groups = "keep") |> 
    dplyr::arrange(t_break)
  
  intervals <- cut(grid_smooth, breaks = c(0, break_points$t_break, 1), 
                   include.lowest = TRUE)
  intervals_H <- intervals
  levels(intervals_H) <- c(break_points$H, break_points$H[length(break_points$H)])
  
  intervals_L <- intervals
  levels(intervals_L) <- c(break_points$L, break_points$L[length(break_points$L)])
  
  intervals_sigma <- intervals
  levels(intervals_sigma) <- c(break_points$sigma, 
                               break_points$sigma[length(break_points$sigma)])
  grid_tibble <- tibble::tibble(t = grid_smooth, H = as.numeric(as.character(intervals_H)),
                                L = as.numeric(as.character(intervals_L)),
                                sigma = as.numeric(as.character(intervals_sigma)))
  #calculate the constants needed to smooth bandwidth 
  #here we use the one associated with epanechnikov kernel
  cst_kernel <- 1.5 * (1/(1 + 2 * grid_tibble$H) - 1/(3 + 2 * grid_tibble$H))
  #cst_kernel <- 1/4 * gamma(grid_tibble$H + 1) + 1/4 * gamma(grid_tibble$H + 2) 
  q1 <- grid_tibble$L/factorial(floor(grid_tibble$H)) * sqrt(cst_kernel)
  q2 <- grid_tibble$sigma
  moments <- estimate_variance_curves(data = curves, 
                                      params = params, 
                                      grid_smooth = grid_smooth, 
                                      sigma = grid_tibble$sigma, 
                                      mu0 = unique(params$mu0))
  q3 <- sqrt(moments$varXt)
  
  M_length <- curves |> purrr::map_dbl(~(length(.x$t)))
  M_length_cum <- cumsum(M_length)
  
  vec <- seq(1, max(M_length_cum))
  vec_cut <- cut(vec, breaks = c(1, M_length_cum), labels = FALSE, 
                include.lowest = TRUE)
  #seq of indices
  Mi_vec <- split(vec, vec_cut) |> unname()
  
  T_diff <- lapply(curves, function(x) x$t) |> 
    unlist() |>
    outer(grid_smooth, FUN = "-")
  
  indic <- T_diff |>
    abs() |>
    outer(grid_bandwidth, FUN = "<=") * 1
  
  #dim G x N x H
  wi <- sapply(Mi_vec, function(Mi) {
    indic[Mi,,] |> colSums(na.rm = TRUE) |>
      (\(x) x >= k0)() * 1
  }, simplify = "array") |>
    aperm(c(1, 3, 2))

  #dim G x H
  WN <- wi |> aperm(c(2, 1, 3)) |> colSums(na.rm = TRUE)
  
  #compute max Wm
  #only compute over max of m
  #dim G x N X H
  max_Wm_num <- sapply(Mi_vec, function(Mi) {
    T_diff[Mi, ] |> abs() |> matrixStats::colMins() 
  }, simplify = "array") |> outer(grid_bandwidth, FUN = "/") |>
    epa_kernel()
  
  K_sum <- T_diff |>
    outer(grid_bandwidth, FUN = "/") |>
    epa_kernel()

  Wm_denom <- sapply(Mi_vec, function(Mi) {
    K_sum[Mi, ,] |> colSums(na.rm = TRUE)
  }, simplify = "array") |> aperm(c(1, 3, 2))

  Wm <- max_Wm_num / Wm_denom
  Wm[is.nan(Wm)] <- 0
  
  #each matrix has dim (G x N)
  Ni <- wi / Wm
  Ni[is.nan(Ni)] <- 0
  
  #dim G x H
  WN_inv2 <- 1/WN^2
  
  wiNi <- wi / Ni
  wiNi[is.nan(wiNi)] <- 0
  #dim G x H
  wiNi_sum <- wiNi |> aperm(c(2, 1, 3)) |> colSums(na.rm = TRUE)
  
  Nmu <- 1/(WN_inv2*wiNi_sum)
  #Nmu[is.nan(Nmu)] <- 0
  
  #column recycling - dim G x H for all q terms
  q2_term <- q2^2 / Nmu
  #q2_term[!is.finite(q2_term)] <- 0
  
  q1_term <- sapply(grid_bandwidth, function(h) {
    q1^2 * h^(2*grid_tibble$H)
  })
  #column recycling 
  q3_term <- q3^2  * (1/WN - 1/length(curves))
  
  risk <- q1_term + q2_term + q3_term
  
  min_h_index <- apply(risk, 1, which.min)
  
  grid_tibble |> dplyr::mutate(bandwidth = grid_bandwidth[min_h_index])
}

#input:
#curves: list of vectors
#params: tibble, with estimate H, L, sigma
#output: G x N matrix 
smooth_curves_mean <- function(curves, params, grid_bandwidth,
                               grid_smooth, k0) {
  
  params <- estimate_bandwidth_mean(curves = curves, params = params, 
                                    grid_bandwidth = grid_bandwidth,
                                    grid_smooth = grid_smooth, 
                                    k0 = k0)

  M_length <- curves |> purrr::map_dbl(~(length(.x$t)))
  M_length_cum <- cumsum(M_length)
  
  vec <- seq(1, max(M_length_cum))
  vec_cut <- cut(vec, breaks = c(1, M_length_cum), labels = FALSE, 
                 include.lowest = TRUE)
  Mi_vec <- split(vec, vec_cut) |> unname()
  
  T_diff <- lapply(curves, function(x) x$t) |> 
    unlist() |>
    outer(grid_smooth, FUN = "-")
  
  indic <- T_diff |>
    abs() |>
    (\(x) x <= matrix(params$bandwidth, nrow = sum(M_length), 
                      ncol = length(grid_smooth), byrow = TRUE))() * 1
  
  #dim G x N 
  wi <- sapply(Mi_vec, function(Mi) {
    indic[Mi,] |> colSums(na.rm = TRUE) |>
      (\(x) x >= k0)() * 1
  })
  
  #dim G x 1
  WN <- wi |> rowSums(na.rm = TRUE)
  
  Wm_num <- (T_diff / matrix(params$bandwidth, nrow = sum(M_length), 
                        ncol = length(grid_smooth), byrow = TRUE)) |>
    epa_kernel()
  
  #transposing makes the next calculation easier due to alignment of 
  #column indices - N x G
  Wm_denom <- sapply(Mi_vec, function(Mi) {
    Wm_num[Mi, ] |> colSums(na.rm = TRUE)
  }) |> t()

  #we repeat the matrix to have the same dimensions as the numerator 
  #for easy element by element division
  Wm_denom_rep <- matrix(data = rep(NA, times = sum(M_length) * length(grid_smooth)),
                         nrow = sum(M_length),
                         ncol = length(grid_smooth))
  
  for(i in seq_along(M_length)) {
    Wm_denom_rep[Mi_vec[[i]], ] <- matrix(Wm_denom[i, ], 
                                   nrow = M_length[i], 
                                   ncol = length(grid_smooth),
                                   byrow = TRUE)
  }
  #dim M x G
  Wm <- Wm_num / Wm_denom_rep
  Wm[is.nan(Wm)] <- 0
  
  Wm_reshaped <- lapply(Mi_vec, function(Mi) {
    Wm[Mi, ]
  })
  #dim G x N
  Xt <- sapply(seq_along(curves), function(i) {
    t(Wm_reshaped[[i]]) %*% curves[[i]]$x
  })
  
  #final output: N x G matrix
  list(params = params, x_hat = t(Xt))
}

#input:
#1. observed curves
#2. tibble of parameters
mean_ll <- function(data, params, grid_bandwidth = lseq(0.005, 
                                                        0.1, length.out = 151),
                    grid_smooth = seq(0, 1, length.out = 101),
                    k0 = 2) {
  smoothed <- smooth_curves_mean(curves = data, params = params, 
                                 grid_bandwidth = grid_bandwidth,
                                 grid_smooth = grid_smooth, k0 = k0)
  mu_hat <- smoothed$x_hat |> colMeans(na.rm = TRUE)
  list(params = smoothed$params, mu_hat = mu_hat)
}


smooth_curves_mean_plugin <- function(curves, params, grid_smooth,
                                      k0, bandwidth) {
  
  M_length <- curves |> purrr::map_dbl(~(length(.x$t)))
  
  #dim Mi x G
  Wm_num <- purrr::map(curves, ~outer(.x$t, grid_smooth, "-")) |>
    lapply(function(Ti) {
      sapply(seq_along(grid_smooth), function(t) {
        outer(Ti[, t], bandwidth[, t], FUN = "/")
      }, simplify = "array") |> epa_kernel()
    })
  
  Wm_denom <- purrr::map(Wm_num, ~colSums(.x, na.rm = TRUE)) |>
    purrr::map2(M_length, ~replicate(.y, .x)) |> 
    purrr::map(~aperm(.x, c(3, 1, 2)))
  
  Wm <- purrr::map2(Wm_num, Wm_denom, ~(.x / .y) |> 
                      (\(x) replace(x, is.nan(x), 0))()) 
  
  Xt <- purrr::map2(curves, Wm, function(i, w) {
    sapply(seq_along(grid_smooth), function(s) {
      t(w[,,s]) %*% i$x
    })
  })
  
  return(Xt)
}
