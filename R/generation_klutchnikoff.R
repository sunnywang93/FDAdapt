## =====================================================
## Generation of MFBM
## See the following paper for the covariance structure
## https://doi.org/10.3390/fractalfract6020074
## Equations (6) and (7)
## =====================================================


##
## Usual hurst functions
##

hurst_atan <- function(t) {
  atan(t) / pi + 0.5
}

hurst_linear <- function(t, h_left = 0.2, h_right = 0.8) {
  t1 <- max(t)
  t0 <- min(t)
  a <- (h_right - h_left) / (t1 - t0)
  b <- h_right - a * t1
  pmin(a * t + b, 1)
}

hurst_logistic <- function(t,
                           h_left = 0.2,
                           h_right = 0.8,
                           slope = 30,
                           change_point_position = 0.5) {
  change_point <- change_point_position * (max(t) + min(t))
  u <- (t - change_point) / (max(t) - min(t))
  (h_right - h_left) / (1 + exp(- slope * u)) + h_left
}


# not totally in agreement with the theory
hurst_piecewise <- function(t,
                            h_left = 0.2,
                            h_right = 0.8,
                            change_point_position = 0.5) {
  change_point <- change_point_position * (max(t) + min(t))
  h_left * (t <= change_point) + h_right * (t > change_point)
}


##
## Covariance structure
##

constant_d <- function(x, y) {
  a <- gamma(2 * x + 1) * gamma(2 * y + 1) * sin(pi * x) * sin(pi * y)
  b <- 2 * gamma(x + y + 1) * sin(pi * (x + y) / 2)
  sqrt(a) / b
}

covariance_mfbm <- function(points, hurst,...) {
  tmp <- expand.grid(s = points, t = points)
  s <- tmp$s
  t <- tmp$t
  hs <- hurst(s, ...)
  ht <- hurst(t, ...)
  hh <- hs + ht
  values <- constant_d(hs, ht) *
    (t**hh + s**hh - abs(t - s)**hh)
  matrix(values, ncol = length(points))
}

covariance_OU <- function(points, l) {
  d <- abs(outer(points, points, FUN = "/"))
  exp(-d/l)
}

##
## data generating process
##

#output: list of time points, one vector for each curve
generates_points <- function(N, m, distribution = runif, ...) {
  M <- rpois(N, m)
  lapply(M, function(x) {
      sort(distribution(x, ...))
    })
}

generate_curves <- function(points_list,
                            hurst,
                            sigma = 0.1,
                            add_regular_grid_of_size = 0,
                            tau,
                            L = 1) {
  # delta <- (t_max - t_min) / 10
  # t_min <- min(sapply(points_list, min)) - delta
  # t_max <- max(sapply(points_list, max)) + delta
  t_min <- 0
  t_max <- 1
  grid <- seq(t_min, t_max, length.out = add_regular_grid_of_size)
  #generate covariance function for each curve
  covariance_list <- lapply(points_list, function(point) {
      pp <- sort(c(point, grid))
      covariance_mfbm(pp, hurst) 
    })


  #generate curves from multivariate normal with covariance function above
  curves <- lapply(covariance_list, function(covariance) {
      MASS::mvrnorm(1, mu = rep(tau * rnorm(1), ncol(covariance)), 
                    Sigma = L * covariance) 
    })
  x_min <- min(sapply(curves, min))
  x_max <- max(sapply(curves, max))
  #extract the curves on the grid and the ones with Mi time points
  curves_random <- list()
  curves_grid <- list()
  for (i in 1:length(curves)) {
    pp <- sort(c(points_list[[i]], grid))
    idx <- which(pp %in% grid)
    curves_grid[[i]] <- curves[[i]][idx]
    curves_random[[i]] <- curves[[i]][-idx]
  }
  out_curves <- list()
  for (i in seq_along(points_list)) {
    ideal <- list(
      t = grid,
      x = curves_grid[[i]]
    )
    observed <- list(
      t = points_list[[i]],
      x = curves_random[[i]] +
        rnorm(length(curves_random[[i]]), 0, sigma)
    )
    class(ideal) <- "curves"
    class(observed) <- "curves"
    out_curves[[i]] <- list(
      ideal = ideal,
      observed = observed
    )
  }
  out <- list()
  out$t_min <- t_min
  out$t_max <- t_max
  out$x_min <- x_min
  out$x_max <- x_max
  out$curves <- out_curves
  class(out) <- "curves_ideal_observed"
  out
}

##
## Class "curves_ideal_observed": methods for generic functions
##

plot.curves_ideal_observed <- function(curves, which_curves = 1:length(curves$curves),
  ...
) {
  courbes <- curves$curves[which_curves]
  my_colors <- rainbow(length(courbes))
  plot(NULL,
       type = "n",
       xlim = c(curves$t_min, curves$t_max),
       ylim = c(curves$x_min, curves$x_max),
       ...
  )
  for (i in seq_along(courbes)) {
    courbe <- courbes[[i]]
    mc <- my_colors[i]
    with(
      courbe,
      {
        points(observed$t, observed$x, col = mc)
        lines(ideal$t, ideal$x, col = mc)
      }
    )
  }
}


plot_subset <- function(curves, t_min, t_max, x_min, x_max, ...) 
  {
  courbes <- curves
  my_colors <- rainbow(length(courbes))
  plot(NULL,
       type = "n",
       xlim = c(t_min, t_max),
       ylim = c(x_min, x_max),
       ...
  )
  for (i in seq_along(courbes)) {
    courbe <- courbes[[i]]
    mc <- my_colors[i]
    with(
      courbe,
      {
        points(observed$t, observed$x, col = mc)
        lines(ideal$t, ideal$x, col = mc)
      }
    )
  }
}



# hurst_path_logistic <- function(t, h_left, h_right, change_point, slope) {
#   hurst_logistic(
#     t,
#     h_left = h_left,
#     h_right = h_right,
#     change_point_position = change_point,
#     slope = slope
#   )
# }
# 
# hurst_path_piecewise <- function(t, h_left, h_right, change_point) {
#   hurst_piecewise(
#     t,
#     h_left = h_left,
#     h_right = h_right,
#     change_point_position = change_point
#   )
# }
# 



