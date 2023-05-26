###############################################################################
### Data generation process
###
### Generation of MFBM
### See the following paper for the covariance structure
### https://doi.org/10.3390/fractalfract6020074
### Equations (6) and (7)
###############################################################################


hurst_atan <- function(t_vec) {
    atan(t_vec) / pi + 0.5
}

#' Generate Hurst function with linear change points
#'
#' Generate linear Hurst functions with two Hurst indexes.
#'
#' @param t_vec Vector of sampling points.
#' @param h_left First Hurst Index.
#' @param h_right Second Hurst Index
#' @returns A vector of length `t_vec` containing the Hurst function.
hurst_linear <- function(t_vec, h_left = 0.2, h_right = 0.8) {
    t1 <- max(t_vec)
    t0 <- min(t_vec)
    a <- (h_right - h_left) / (t1 - t0)
    b <- h_right - a * t1
    pmin(a * t_vec + b, 1)
}


#' Generate Hurst function with logistic change points
#'
#' Generate logistic Hurst functions with two Hurst indexes.
#'
#' @param t_vec Vector of sampling points.
#' @param h_left First Hurst Index.
#' @param h_right Second Hurst Index
#' @param slope Slope of transition between `h_left` and `h_right`.
#' @param change_point_position Change point position between two Hurst indexes.
#' @returns A vector of length `t_vec` containing the Hurst function.
#' @export
hurst_logistic <- function(
    t_vec,
    h_left = 0.5,
    h_right = 0.8,
    slope = 15,
    change_point_position = 0.5
) {
    change_point <- change_point_position * (max(t_vec) + min(t_vec))
    u <- (t_vec - change_point) / (max(t_vec) - min(t_vec))
    (h_right - h_left) / (1 + exp(-slope * u)) + h_left
}

#' Generate piecewise Hurst function.
#'
#' Generate piecewise Hurst functions with two Hurst indexes. Will have
#' discontinuous jumps, use with caution.
#'
#' @param t_vec Vector of sampling points.
#' @param h_left First Hurst Index.
#' @param h_right Second Hurst Index
#' @param change_point_position Change point position between two Hurst indexes.
#' @returns A vector of length `t_vec` containing the Hurst function.
hurst_piecewise <- function(
    t_vec,
    h_left = 0.2,
    h_right = 0.8,
    change_point_position = 0.5
) {
    change_point <- change_point_position * (max(t_vec) + min(t_vec))
    h_left * (t_vec <= change_point) + h_right * (t_vec > change_point)
}


#' Calculate constants associated to multi-fractional brownian motion.
#'
#' @param x Vector of Hurst index in one dimension.
#' @param y Vector of Hurst index in second dimension.
#' @returns A vector of length `x`/ `y` containing the constants.
#' @export
constant_d <- function(x, y) {
    a <- gamma(2 * x + 1) * gamma(2 * y + 1) * sin(pi * x) * sin(pi * y)
    b <- 2 * gamma(x + y + 1) * sin(pi * (x + y) / 2)
    sqrt(a) / b
}

#' Generate covariance function associated to multi-fractional brownian motion.
#'
#' @param points Vector of sampling points.
#' @param hurst Vector containing the Hurst function.
#' @param points_disto Vector containing the distorted time points.
#' @param norm TRUE/FALSE Whether to normalise the covariance by constants.
#' @returns A vector of length `x`/ `y` containing the constants.
#' @export
covariance_mfbm <- function(points, hurst, points_disto, ..., norm = FALSE) {
    tmp <- expand.grid(s = points, t = points)
    s <- tmp$s
    t <- tmp$t
    hs <- hurst(s)
    ht <- hurst(t)
    hh <- hs + ht

    tmp <- expand.grid(s = points_disto, t = points_disto)
    s_ <- tmp$s
    t_ <- tmp$t
    values <- constant_d(hs, ht) * (t_**hh + s_**hh - abs(t_ - s_)**hh)
    if (norm) {
        values <- values / sqrt(s_**hh * t_**hh)
    }
    values[!is.finite(values)] <- 0
    matrix(values, ncol = length(points))
}

##
##  Estimate derivatives
##
estimate_derivative <- function(crv, #beta, const,
                                grid_size = 101L, range.x = c(0,1)) {
    curves <- crv
    for (i in seq_along(curves)) {
        t <- curves[[i]]$t
        x <- curves[[i]]$x
        m <- length(t)
        h <- 3 / m #log(m) / m #(1/m)**(1/(2*beta+1))*const
        LP <- KernSmooth::locpoly(
            t,
            x,
            drv = 1,
            degree = 2,
            bandwidth = h,
            gridsize = grid_size,
            range.x = range.x
        )
        names(LP) <- c("t", "x")
        curves[[i]] <- LP
    }
    #structure(curves, class = "curves_with_estimated_derivatives")
    curves
}


#' Generate random time points for functional data.
#'
#' Generates sampling points in the random design framework, where the exact
#' number of points per curve are generated according to a Poisson distribution
#' with `m` average number of points.
#'
#' @param N Number of curves.
#' @param m Average number of observed points per curve.
#' @param distribution The distribution of the sampling points. Defaults to the
#' uniform distribution.
#' @returns A list containing the generated time points for `N` curves.
#' @export
generate_points <- function(N, m, distribution = runif, ...) {
    M <- rpois(N, m)
    lapply(M, function(x) {
            sort(distribution(x, ...))
      })
}

#' Generate curves in the functional data framework.
#'
#' @param points_list List where each element contains the sampling points
#' for one curve.
#' @param hurst Hurst function.
#' @param distortion_model Distortion function `A(.)` for the time points.
#' @param variance_fun Variance function `v(t)^2`.
#' @param sigma0 Numeric, the baseline noise level of curves.
#' @param hetero Boolean, indicating whether to generate heteroscedastic curves.
#' See `sigma_het` function for more details.
#' @param regular_grid Vector, containing the grid on which the true curves
#' will lie.
#' @param add_one_to_hurst Boolean, for use if fractional regularity (1 + ...)
#' is desired.
#' @param norm_cov Boolean, indicating whether covariance should be normalised
#' by constant associated to multi-fractional brownian motion.
#' @returns A list, containing the ideal and sampled curves, together with other
#' auxiliary information.
#' @export
generate_curves <- function(
    points_list,
    hurst,
    distortion_model = function(x) x,
    variance_fun = function(x) 1,
    sigma0 = 0.1,
    hetero = TRUE,
    regular_grid = seq(0, 1, l = 101),
    add_one_to_hurst = FALSE,
    norm_cov = FALSE,
    L = 1,
    novar = FALSE,
    ...
) {
    grid <- regular_grid
    # Create covariance function with parameter values
    covariance_list <- lapply(points_list, function(point) {
            pp <- sort(c(point, grid))
            pp_disto <- distortion_model(pp)
            covariance_mfbm(pp, hurst, pp_disto, norm = norm_cov)
            })

    # Simulate values from covariance function
    curves <- lapply(covariance_list, function(covariance) {
      MASS::mvrnorm(1, mu = rep(0, ncol(covariance)), Sigma = L * covariance)
      })

    # Create new sample paths with matched variance
    if(novar == FALSE) {
      curves <- lapply(1:length(curves), function(idx) {
        pp <- sort(c(points_list[[idx]], grid))
        sqrt(variance_fun(pp)) * distortion_model(pp)**(-hurst(pp)) * curves[[idx]]
      })
    }
    # Add random starting point for curves
    # curves <- lapply(curves, function(curve) {
    #   curve + tau * rnorm(1)
    #   })

    x_min <- min(sapply(curves, min))
    x_max <- max(sapply(curves, max))

    curves_random <- list()
    curves_grid <- list()
    curves_deg_0 <- list()
    for (i in 1:length(curves)) {
        pp <- sort(c(points_list[[i]], grid))
        idx <- which(pp %in% grid)
        if (add_one_to_hurst) {
          dt <- diff(c(0, pp))
          curves_deg_0[[i]] <- curves[[i]]
          curves[[i]] <- cumsum(curves_deg_0[[i]] * dt)
          curves_deg_0[[i]] <- curves_deg_0[[i]][idx]
        }
        curves_grid[[i]] <- curves[[i]][idx]
        curves_random[[i]] <- curves[[i]][-idx]
    }

    out_curves <- list()
    for (i in seq_along(points_list)) {
        ideal <- list(
            t = grid,
            x = curves_grid[[i]]
        )
        class(ideal) <- "curves"
        observed <- list(
            t = points_list[[i]],
            x = curves_random[[i]] +
                rnorm(n = length(points_list[[i]]),
                      sd = sigma_het(sigma0, hetero, points_list[[i]]))
        )
        class(observed) <- "curves"
        if (add_one_to_hurst) {
          deg_0 <- list(
            t = grid,
            x = curves_deg_0[[i]]
          )
          class(deg_0) <- "curves"
        } else {
          deg_0 = NULL
        }
        out_curves[[i]] <- list(
            deg_0 = deg_0,
            ideal = ideal,
            observed = observed
        )
    }

    out <- list()
    out$t_min <- min(grid)
    out$t_max <- max(grid)
    out$x_min <- x_min
    out$x_max <- x_max
    out$curves <- out_curves
    class(out) <- "curves_ideal_observed"
    out
}


#' Generate a mean curve.
#'
#' Generates a mean curve using a combination of sin basis functions.
#'
#' @param grid Vector of sampling points to generate the mean curve.
#' @param K Number of basis functions to use.
#' @param alpha Parameter controlling the regularity of the mean curve. A higher
#' `alpha` will produce smoother mean curves.
#' @returns A list, containing the sampling points and observed points.
#' @export

generate_mean_curve <- function(grid = seq(0, 1, length.out = 101), K = 60,
    alpha = 1) {

    K_seq <- 1:K
    Z <- rnorm(K)
    xi <- 1/((K_seq - 1/2)**alpha * pi**2)

    mu_kt <- sapply(grid, function(t) {
        sqrt(2) * Z * sqrt(xi) * sin((K_seq - 0.5) * pi * t)
    })

    mu_t <- colSums(mu_kt)
    list(t = grid, x = mu_t)
}



sample_curve <- function(curve, grid) {
    idx <- map_dbl(grid, ~which.min(abs(.x - curve$t)))
    curve$x[idx]
}


#' Generates gaussian noise for function data
#'
#' Generates possibly heteroscedastic gaussian noise of functional data,
#' based on the observed time points of one curve.
#' Heteroscedasticity is generated using sin functions.
#'
#' @param sigma0 Numeric, baseline noise level where the noise will be centered
#' around.
#' @param hetero Boolean, where TRUE indicate
#' @param t_vec Vector, containing the sampling points for
#' @returns A list containing the noise for each curve if `hetero = TRUE`. If
#' not, a numeric will be returned.
#' @export
sigma_het <- function(sigma0, hetero = TRUE, t_vec) {
  if(hetero) {
    sigma0 * (1 + sin(8 * pi * t_vec) / 3)
  } else {
    sigma0
  }
}




#' Plot functions for generated true curves
#'
#' @param curves List, containing the ideal curves.
#' @param which_curves Vector, indicating the curves to be plotted.
#' @export
plot.curves_ideal_observed <- function(
        curves,
        which_curves = 1:length(curves$curves),
        ...
) {
    courbes <- curves$curves[which_curves]
    my_colors <- rainbow(length(courbes), alpha = 0.5)
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

plot_derivatives <- function(
    curves,
    which_curves = 1:length(curves$curves),
    ...
) {
  courbes <- curves$curves[which_curves]
  my_colors <- rainbow(length(courbes), alpha = 0.5)
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
        lines(deg_0$t, deg_0$x, col = mc)
        lines(estimated_derivative$t, estimated_derivative$x, col = mc)
      }
    )
  }
}

plot.curves_with_estimated_derivatives <- function(
        curves,
        which_curves = 1:length(curves$curves),
        ...
) {
    courbes <- curves$curves[which_curves]
    my_colors <- rainbow(length(courbes), alpha = 0.5)
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
                lines(deg_0$t, deg_0$x, col = mc)
                lines(estimated_derivative$t, estimated_derivative$x, col = mc)
            }
        )
    }
}

