
library(MASS)
library(rpart)
library(matrixStats)
library(tibble)
library(dplyr)

#sigma estimation is bad - can't seem to figure out why
estimate_sigma <- function(data) {
  Mbar <- data |> purrr::map_dbl(~length(.x$t)) |> mean() 
  delta <- (1/log(Mbar)^2)
  indic <- data |> 
    purrr::map(~as.numeric(abs(diff(sort(.x$t, decreasing = TRUE))) <= delta))
  diffsq <- data |> purrr::map(~diff(sort(.x$x, decreasing = TRUE))^2)
  sum_diffsq <- purrr::map2_dbl(diffsq, indic, ~sum(.x * .y))
  denom <- indic |> purrr::map_dbl(~sum(.x))
  sum_diffsq_norm <- sum_diffsq / (2 * denom)
  sqrt(mean(sum_diffsq_norm))
}



estimate_mu0 <- function(data) {
  T_all <- data %>% purrr::map(~.x$t) %>% unlist() %>% sort()
  min(density(T_all, from = 0.1, to = 0.9)$y)
}




presmoothing_piecewise <- function (data, t2 = seq(.2, .8, l = 20), 
                                    drv = 0, init_b, sigma = NULL, L, mu0 = NULL)
{
  NW <- function(t, T_mi, Y_mi, h, alpha = init_b) {
    K <- function(x, beta){
      ((1 + beta)/ (2 * beta))  * (1 - abs(x)^beta) * (abs(x) <= 1)
    }
    tt <- matrix(rep(t, length(T_mi)), ncol = length(T_mi))
    TT <- t(matrix(rep(T_mi, length(t)), ncol = length(t)))
    Kx <- 2 * K(x = (tt - TT) / h, beta = alpha)
    (Kx %*% Y_mi) / matrix(rowSums(Kx))
  }
  m <- data %>% purrr::map_dbl(~length(.x$t)) %>% mean()
  delta <- min(log(m)^(-1.1), 0.2)
  t1 <- t2 - delta/2
  t3 <- t2 + delta/2
  t_list <- rbind(t1, t2, t3)
  if(is.null(sigma)) {
    sigma <- estimate_sigma(data)
  }
  if(is.null(mu0)) {
    mu0 <- estimate_mu0(data)
  }
  #init_b and L must be the same length as t2
  if(length(init_b) == 1) {
    init_b <- rep(init_b, length(t2))
  } 
  if(length(L) == 1) {
    L <- rep(L, length(t2))
  }
  
  c <- (sigma^(2*init_b) * L * ((init_b + 1)/ 2 * init_b^2 * mu0)^init_b)^(1 / (2*init_b + 1))
  psi_m <- (1 / m)^(init_b / (2 * init_b + 1))
  b_naive <- pmax(pmin((c * psi_m / L)^(1 / init_b), delta/4), log(m)/m)

  inner_loop <- function(i, data, t_list, b_naive, init_b) {
    sapply(data, function(x) {
      NW(t = t_list[, i],
         T_mi = x$t, Y_mi = x$x,
         h = b_naive[i], alpha = init_b[i])
    }) %>% t() 
  }
  
  purrr::map(1:ncol(t_list), ~list(t_list = t_list[,.x],
                                   x = inner_loop(.x, data = data, t_list = t_list,
                                                  b_naive = b_naive,
                                                  init_b = init_b)))
}


estimate_theta <- function(data, center = TRUE) {
  data %>% purrr::map_dbl(function(d) {
    if(center) {
      a <- mean((d$x[, 1] - d$x[, 3] - mean(d$x[, 1] - d$x[, 3], na.rm = TRUE))^2, na.rm = TRUE)
      b1 <- mean((d$x[, 2] - d$x[, 3] - mean(d$x[, 2] - d$x[, 3], na.rm = TRUE))^2, na.rm = TRUE)
      b2 <- mean((d$x[, 1] - d$x[, 2] - mean(d$x[, 1] - d$x[, 2], na.rm = TRUE))^2, na.rm = TRUE)
    }
    else {
      a <- mean((d$x[, 1] - d$x[, 3])^2, na.rm = TRUE)
      b1 <- mean((d$x[, 2] - d$x[, 3])^2, na.rm = TRUE)
      b2 <- mean((d$x[, 1] - d$x[, 2])^2, na.rm = TRUE)
    }
    min(max(((log(a) - log(b1)) + (log(a) - log(b2))) / (4*log(2)), 0.1, na.rm = TRUE),
        1, na.rm = TRUE)
  })
}

estimate_L <- function(data, H0_list) {
  L1 <- data %>% purrr::map2(H0_list, ~(mean((.x$x[, 2] - .x$x[, 3] - mean(.x$x[, 2] - .x$x[, 3], na.rm = TRUE))^2,
                                             na.rm = TRUE)/abs(.x$t[3] - .x$t[2])^(2*.y)))
  L2 <- data %>% purrr::map2(H0_list, ~(mean((.x$x[, 1] - .x$x[, 2] - mean(.x$x[, 1] - .x$x[, 2], na.rm = TRUE))^2,
                                             na.rm = TRUE)/abs(.x$t[2] - .x$t[1])^(2*.y)))
  L1_vec <- unlist(L1, use.names = FALSE)
  L1_vec[is.nan(L1_vec)] <- 0
  L2_vec <- unlist(L2, use.names = FALSE)
  L2_vec[is.nan(L2_vec)] <- 0
  as.vector(sqrt((L1_vec + L2_vec)/2))
}


L_tree <- function(L, H_tree) {
  tibble(H = H_tree, L = L) %>% group_by(H) %>% mutate(L_tree = median(L)) %>%
    ungroup() %>% pull(L_tree)
}

estimate_holder_const <- function(data, sigma = NULL, mu0 = NULL,
                                  grid_estim = seq(0.2, 0.8, length.out = 20)) {
  if(is.null(sigma)) {
    sigma <- estimate_sigma(data)
  }
  if(is.null(mu0)) {
    mu0 <- estimate_mu0(data)
  }
  presmoothed <- data %>% presmoothing_piecewise(init_b = 1, L = 1, sigma = sigma, mu0 = mu0,
                                                 t2 = grid_estim)
  H0 <- estimate_theta(presmoothed)
  L0 <- estimate_L(presmoothed, H0)
  tree1_tib <- tibble(y = H0, x = grid_estim)
  tree1 <- partykit::as.party(rpart(y ~ x, data = tree1_tib))
  H0_tree <- as.vector(predict(tree1, x = grid_estim, FUN = function(y, w) median(y)))
  L0_tree <- L_tree(L0, H0_tree)
  
  presmoothed_2nd <- data %>% presmoothing_piecewise(init_b = H0_tree, L = L0_tree, sigma = sigma, mu0 = mu0,
                                                     t2 = grid_estim)
  H1 <- estimate_theta(presmoothed_2nd)
  L1 <- estimate_L(presmoothed_2nd, H1)
  tree2_tib <- tibble(y = H1, x = grid_estim)
  tree2 <- partykit::as.party(rpart(y ~ x, data = tree2_tib))
  H1_tree <- as.vector(predict(tree2, x = grid_estim, FUN = function(y, w) median(y)))
  L1_tree <- L_tree(L1, H1_tree)
  tibble(t = grid_estim, H = H1_tree, L = L1_tree, sigma = sigma, mu0 = mu0)
}


