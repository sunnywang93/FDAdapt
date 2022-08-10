source("/Users/swang/Dropbox/research_code/mean_optimised.R")

estimate_bandwidth_covariance <- function(curves, params, grid_bandwidth,
                          grid_smooth, k0)
{
  break_points <- params |> dplyr::group_by(H, L, sigma) |>
    dplyr::summarise(t_break = max(t), .groups = "keep") |> 
    dplyr::arrange(t_break)
  
  
  intervals <- cut(grid_smooth, breaks = c(0, break_points$t_break, 1), 
                   include.lowest = TRUE)
  intervals_H <- intervals
  levels(intervals_H) <- c(break_points$H, 
                           break_points$H[length(break_points$H)])
  
  intervals_L <- intervals
  levels(intervals_L) <- c(break_points$L, break_points$L[length(break_points$L)])
  
  intervals_sigma <- intervals
  levels(intervals_sigma) <- c(break_points$sigma, 
                               break_points$sigma[length(break_points$sigma)])
  grid_tibble <- tibble::tibble(t = grid_smooth, H = as.numeric(as.character(intervals_H)),
                        L = as.numeric(as.character(intervals_L)),
                        sigma = as.numeric(as.character(intervals_sigma)))
  
  cst_kernel <- 1/4 * gamma(grid_tibble$H + 1) + 1/4 * gamma(grid_tibble$H + 2) 
  var_ <- estimate_variance_curves(data = curves, params = params, 
                                   grid_smooth = grid_smooth,
                                   sigma = unique(params$sigma), 
                                   mu0 = unique(params$mu0))
  # q1 <- sqrt(2 * var_$EXt2) * sqrt(cst_kernel) / factorial(floor(grid_tibble$H)) *
  #   grid_tibble$L
  q1 <- sqrt(2 * var_$EXt2 * cst_kernel / factorial(floor(grid_tibble$H))^2 *
                 grid_tibble$L)
  
  q2 <- grid_tibble$sigma * sqrt(var_$EXt2)
  
  q3 <- sqrt(var_$varXtXs / 2)
  
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
  ws <- sapply(Mi_vec, function(Mi) {
    indic[Mi,,] |> colSums(na.rm = TRUE) |>
      (\(x) x >= k0)() * 1
  }, simplify = "array") |>
    aperm(c(1, 3, 2))
  
  #dim G x G x N x H
  wst <- apply(ws, MARGIN = c(2, 3), function(x) x %*% t(x)) |>
    array(dim = c(length(grid_smooth), 
                  length(grid_smooth), 
                  length(curves),
                  length(grid_bandwidth)))
  

  #dim G x G x H
  WN <- wst |> aperm(c(3, 4, 1, 2)) |> 
    colSums(na.rm = TRUE) |> 
    aperm(c(2, 3, 1))
  
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
  
  rm(K_sum)
  
  Wm <- max_Wm_num / Wm_denom
  Wm[is.nan(Wm)] <- 0
  
  #replicate into new 4th dimension s
  Wm_rep <- replicate(length(grid_smooth), Wm) |>
    aperm(c(4,1,2,3))
  
  Nt <- wst / Wm_rep
  
  rm(Wm_rep)
  
  Nt[is.nan(Nt)] <- 0
  
  gamma_sum <- (wst / Nt) |> aperm(c(3,1,2,4)) |> colSums(na.rm = TRUE) 
  
  rm(wst, Nt)
  
  Ngamma <- 1 / ((1/WN^2) * gamma_sum)
  #column recycling - q2(s|t)
  q2_term <- array(q2^2, dim = c(length(grid_smooth), 
                                 length(grid_smooth), 
                                 length(grid_bandwidth))) |> 
    aperm(c(2,1,3)) |>
    (\(x) x / Ngamma)()
  #q1(s|t)
  q1_term <- sapply(grid_bandwidth, function(h) {
    q1^2 * h^(2*grid_tibble$H)
  }) |> 
    (\(x) replicate(length(grid_smooth), x))() |>
    aperm(c(1, 3, 2))
    
  q3_term <- replicate(length(grid_bandwidth), q3^2) * (1/WN - 1/length(curves))
  
  risk_s <- q1_term + q2_term + q3_term
  #risk(t|s) is the transpose of risk(s|t)
  risk <- risk_s + aperm(risk_s, c(2, 1, 3))
  
  min_h_index <- apply(risk, MARGIN = c(1,2), which.min)
  
  apply(min_h_index, 2, function(id) grid_bandwidth[id])
}

smooth_curves_covariance <- function(curves, params, grid_bandwidth = 
                                       lseq(0.005, 0.1, length.out = 151),
                                     grid_smooth = seq(0, 1, length.out = 101),
                                     k0 = 2) {
  
  
  bandwidth_matrix <- estimate_bandwidth_covariance(curves, params,
                                                    grid_bandwidth,
                                                    grid_smooth,
                                                    k0)
  
  M_length <- curves |> purrr::map_dbl(~(length(.x$t)))
  
  #dim Mi x G
  Wm_num <- purrr::map(curves, ~outer(.x$t, grid_smooth, "-")) |>
    lapply(function(Ti) {
      sapply(seq_along(grid_smooth), function(t) {
        outer(Ti[, t], bandwidth_matrix[, t], FUN = "/")
      }, simplify = "array") |> epa_kernel()
    })

  Wm_denom <- purrr::map(Wm_num, ~colSums(.x, na.rm = TRUE)) |>
    purrr::map2(M_length, ~replicate(.y, .x)) |> 
    purrr::map(~aperm(.x, c(3, 1, 2)))

  Wm <- purrr::map2(Wm_num, Wm_denom, ~(.x / .y) |> 
                      (\(x) replace(x, is.nan(x), 0))()) 

  XtXs <- purrr::map2(curves, Wm, function(i, w) {
    sapply(seq_along(grid_smooth), function(s) {
      t(w[,,s]) %*% i$x
    })
  }) |> purrr::map(~(.x * t(.x)))
  
    
  list(prod = XtXs, bw = bandwidth_matrix)
}

#center: choose between
# 1. curves (center at the beginning)
# 2. prod_hmu (center with prod (mu bandwidth))
# 3. prod_hgamma (center with prod (gamma bandwidth))

covariance_ll <- function(curves, params, grid_bandwidth = 
                          lseq(0.005, 0.1, length.out = 151),
                          grid_smooth = seq(0, 1, length.out = 101),
                          k0 = 2, center = "curves") {
  
  if(center == "curves") {
    mu <- mean_ll(curves, params, grid_bandwidth, grid_smooth, k0)
    idx <- lapply(curves, function(i) {
      purrr::map_dbl(i$t, ~which.min(abs(.x - mu$params$t)))
    })
    mu_Ti <- purrr::map(idx, ~mu$mu_hat[.x])
    curves <- purrr::map2(curves, mu_Ti, 
                          function(x, mu) {x$x <- x$x - mu; x; })
  }
  
  
  prod_ <- smooth_curves_covariance(curves, params, grid_bandwidth,
                                    grid_smooth, k0)

  wt_cond <- purrr::map(curves, ~abs(outer(.x$t, grid_smooth, "-"))) |>
    lapply(function(Ti) {
      sapply(seq_along(grid_smooth), function(t) {
        (outer(Ti[, t], prod_$bw[, t], FUN = "<=") * 1) |> colSums(na.rm = TRUE)
      }) |>
        (\(x) (x >= k0)  * 1)() 
    }) 
  
  wtws <- purrr::map(wt_cond, ~(.x * t(.x)))

  WN <- Reduce('+', wtws)
  
  prod_sum <- purrr::map2(wtws, prod_$prod, ~(.x * .y)) |>
    (\(x) Reduce('+', x))()

  gamma <- (1/WN) * prod_sum
  
  if(center == "prod_hmu") {
    mu <- mean_ll(curves, params, grid_bandwidth, grid_smooth, k0)
    Gamma <- gamma - (mu$mu_hat %*% t(mu$mu_hat))
  }
  else if(center == "prod_hgamma") {
    Xt_cond <- smooth_curves_mean_plugin(curves, params, grid_smooth,
                                    k0, prod_$bw)
    mu_ts <- purrr::map2(Xt_cond, wt_cond, ~.x * .y) |>
      (\(x) Reduce('+', x) / WN)() |>
      (\(x) x * t(x))()
    Gamma <- gamma - mu_ts
  }
  else {
    Gamma <- gamma
  }
  

  #final step: replace diagonal band
  for (t in 1:ncol(Gamma)) {
    s <- 1
    current_cov <- Gamma[s, t - s + 1]
    while (s <= (t - s + 1)) {
      if (abs(grid_smooth[s] - grid_smooth[t - s + 1]) > 
          prod_$bw[s, t - s + 1]) {
        current_cov <- Gamma[s, t - s + 1]
      } else {
        Gamma[s, t - s + 1] <- current_cov
      }
      s <- s + 1
    }
  }
  #loop over lower anti-diagonal 
  for (s in 1:nrow(Gamma)) {
    t <- ncol(Gamma)
    current_cov <- Gamma[ncol(Gamma) + s - t, t]
    while (t >= (ncol(Gamma) + s - t)) {
      if (abs(grid_smooth[ncol(Gamma) + s - t] - grid_smooth[t]) >
          prod_$bw[ncol(Gamma) + s - t, t]) {
        current_cov <- Gamma[ncol(Gamma) + s - t, t]
      } else {
        Gamma[ncol(Gamma) + s - t, t] <- current_cov
      }
      t <- t - 1
    }
  }
    list(Gamma = Gamma, bw = prod_$bw)
}

