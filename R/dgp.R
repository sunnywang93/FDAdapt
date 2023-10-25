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


#' Generate sampling points for functional data
#'
#' Generates sampling points either for common or random design functional
#' data. For random design, the exact
#' number of points per curve are generated according to a Poisson distribution
#' with `m` average number of points.
#'
#' @param N Number of curves.
#' @param m Average number of observed points per curve.
#' @param distribution The distribution of the sampling points. Defaults to the
#' uniform distribution.
#' @param common Boolean, indicating whether sampling scheme should be common
#' design or sampling design.
#' @param tmin Numeric, left endpoint of the sampling grid.
#' @param tmax Numeric, right endpoint of the sampling grid.
#' @returns A list or a vector, depending on `common`. If `common = FALSE`, a list
#' containing the generated time points for each curve is returned.
#' If `common = TRUE`, a vector containing the sampling points is returned.
#' @export
generate_points <- function(N, m, distribution = runif,
                            common = FALSE, tmin = NULL,
                            tmax = NULL, ...) {

  if(common) {
    if(is.null(tmin) | is.null(tmax)) {
      stop("tmin and tmax must be supplied when common = TRUE")
    } else {
      seq(tmin, tmax, length.out = m)
    }
  } else {
    M <- rpois(N, m)
    lapply(M, function(x) {
      sort(distribution(x, ...))
    })
  }

}

#' Generate random design curves in the functional data framework.
#'
#' @param points_list List where each element contains the sampling points
#' for one curve.
#' @param hurst Hurst function.
#' @param distortion_model Distortion function `A(.)` for the time points.
#' @param variance_fun Variance function.
#' @param sigma0 Numeric, the baseline noise level of curves.
#' @param hetero Boolean, indicating whether to generate heteroscedastic curves.
#' See `sigma_het` function for more details.
#' @param regular_grid Vector, containing the grid on which the true curves
#' will lie.
#' @param add_one_to_hurst Boolean, for use if fractional regularity (1 + ...)
#' is desired.
#' @param norm_cov Boolean, indicating whether covariance should be normalised
#' by constant associated to multi-fractional brownian motion.
#' @param novar Boolean, indicating whether variance function should not be
#' matched. `TRUE` results in variance not being matched.
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
      MASS::mvrnorm(1, mu = rep(0, ncol(covariance)), Sigma = covariance)
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


#' Generate common design curves in the functional data framework.
#'
#' @param points_vector Vector containing the common sampling points for all
#' curves.
#' @param hurst Hurst function.
#' @param distortion_model Distortion function `A(.)` for the time points.
#' @param variance_fun Variance function.
#' @param N Number of curves to simulate.
#' @param sigma0 Numeric, the baseline noise level of curves.
#' @param hetero Boolean, indicating whether to generate heteroscedastic curves.
#' See `sigma_het` function for more details.
#' @param novar Boolean, indicating whether variance function should not be
#' matched. `TRUE` results in variance not being matched.
#' @param norm_cov Boolean, indicating whether covariance should be normalised
#' by constant associated to multi-fractional brownian motion.
#' @returns A list, containing the ideal and sampled curves, together with other
#' auxiliary information.
#' @export
generate_curves_common <- function(points_vector,
                                   hurst,
                                   distortion_model = function(x) x,
                                   variance_fun = function(x) 1,
                                   N,
                                   sigma0 = 0.1,
                                   hetero = TRUE,
                                   novar = FALSE,
                                   norm_cov = FALSE,
                                   ...)
{

  pp_disto <- distortion_model(points_vector)
  covariance <- covariance_mfbm(points_vector,
                                hurst,
                                pp_disto,
                                norm = norm_cov)

  curves <- MASS::mvrnorm(N,
                          mu = rep(0, ncol(covariance)),
                          Sigma = covariance)

  if(novar == FALSE) {
    curves <- apply(curves, 1, function(z) {
      sqrt(variance_fun(points_vector)) *
        distortion_model(points_vector)**(-hurst(points_vector)) *
        z
    })
  }

  out_ideal <- lapply(seq_len(N), function(n) {
    list(t = points_vector,
         x = curves[n, ])
  })

  out_observed <- lapply(seq_len(N), function(n) {
    list(t = points_vector,
         x = curves[n, ] +
           rnorm(n = length(points_vector),
                 sd = sigma_het(sigma0, hetero, points_vector)))
  })


  out <- list()
  out$t_min <- min(points_vector)
  out$t_max <- max(points_vector)
  out$x_min <- min(curves)
  out$x_max <- max(curves)
  class(out) <- "curves_ideal_observed"
  out$curves <- list(ideal = out_ideal,
                     observed = out_observed)
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


#' Generates gaussian standard deviation for function data
#'
#' Generates possibly heteroscedastic standard deviation for gaussian nosie of
#' functional data, based on the observed time points of one curve.
#' Heteroscedasticity is generated using sin functions.
#'
#' @param sigma0 Numeric, baseline noise level where the noise will be centered
#' around.
#' @param hetero Boolean, where TRUE indicate
#' @param t_vec Vector, containing the evaluation points
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


#' Approximate input functions using Fourier bases
#'
#' Performs approximation of input functions (primarily used for eigenfunctions
#' ) using orthonormal Fourier bases and LASSO regression.
#'
#' @param psi_t Vector containing the evaluation points.
#' @param psi_x Vector containing the observed points.
#' @param K Numeric indicating the number of Fourier bases functions to use
#' for approximating input function.
#' @returns Vector containing the coefficients from LASSO regression.
#' @export

psi_approx <- function(psi_t, psi_x, K) {

  # Create orthonormal fourier bases for approximation
  fourier_base <- sqrt(2) * (sin(2 * pi * outer(psi_t, seq_len(K))) +
  cos(2 * pi * outer(psi_t, seq_len(K))))

  # Perform cross-validation to obtain smoothing parameter for LASSO
  pen_cv <- glmnet::cv.glmnet(x = fourier_base,
                      y = psi_x,
                      alpha = 1,
                      type.measure = "mse")

  # Perform LASSO regression with mse-min parameter
  pen_reg <- glmnet::glmnet(x = fourier_base,
                            y = psi_x,
                            alpha = 1,
                            lambda = pen_cv$lambda.min)

  # Bind intercept and coefficients for return
  unname(c(pen_reg$a0, pen_reg$beta[, 1]))

}

# Be careful - K + 1 must be greater than or equal to number of columns of psi!

#' Generates functional data using Karhunen-Loève decomposition
#'
#' Sample paths of a stochastic process are generated using the
#' Karhunen-Loève decomposition, for given eigenvalues and
#' eigenfunction pairs. The input eigenfunctions are approximated
#' by orthogonal Fourier bases using a LASSO regression, before a
#' gram-schmidt orthonormalisation is performed so that the
#' approximated eigenfunctions remain orthonormal.
#'
#' @param lambda Vector of eigenvalues.
#' @param psi List, containing the following elements:
#' - **t** Vector of evaluation points for the eigenfunctions.
#' - **x** Matrix, with each column containing the j-th eigenfunction
#' at evaluation points `t`.
#' @param K Numeric, the number of Fourier bases used to approximate
#' the eigenfunctions. K should be greater than the number of input
#' eigenfunctions you are trying to approximate.
#' @param points_list List, containing the evaluation points for each
#' curve.
#' @param mu List, containing the following elements:
#' - **t** Vector of evaluation points for the mean function.
#' - **x** Vector containing the observed points of the mean function.
#' @param sigma_sd Vector / Numeric, containing the standard deviation
#' structure of the gaussian noise to be added to sample paths.
#' @param xi_fam Character, indicating the distribution to generate
#' the scores from. Gaussian, uniform and chi-squared distributions
#' are currently available.
#' @returns List, with each element containing the following:
#' -**curves** List of curves, containing a vector of evaluation points
#' and observed points.
#' -**mu** List, containing a vector of evaluation points and
#' observed points. The observed points are computed based on the closest point
#' from the mean function on the input evaluation points and the evaluation point
#' of the curves.
#' -**zeta** Vector of length K, containing the true scores.
#' -**v** Matrix with the columns containing the approximated jth eigenfunction.
#' -**sigma** Vector containing the gaussian noise added each curve.
#' @export

generate_KL <- function(lambda,
                        psi,
                        K,
                        points_list,
                        mu,
                        sigma0,
                        hetero,
                        zeta_fam = c("gaussian", "unif", "chi2")) {

  # Approximate each eigenfunction with fourier bases function
  coef_lasso <- apply(psi$x, 2, function(x) psi_approx(psi$t, x, K))

  # Orthonormalise coefficient matrix
  ortho_coef <- pracma::gramSchmidt(A = coef_lasso)

  # Construct basis functions on (possibly) random grid
  fourier_base_grid <- lapply(points_list, function(points) {
    sqrt(2) * (sin(2 * pi * outer(points, seq_len(K))) +
                 cos(2 * pi * outer(points, seq_len(K))))
  }) |>
    purrr::map(~cbind(rep(1, nrow(.x)), .x))

  # Generate scores
  zeta <- lapply(seq_along(points_list), function(i) {
    switch(zeta_fam,
           "gaussian" = sqrt(lambda) * rnorm(n = ncol(psi$x)),
           "unif" = sqrt(lambda) * runif(n = ncol(psi$x),
                                         min = -sqrt(3),
                                         max = sqrt(3)),
           "chi2" = sqrt(lambda / 2) * (rchisq(n = ncol(psi$x),
                                              df = 1) - 1)
           )
  })


  # Simulate sample paths
  paths <- purrr::map2(fourier_base_grid, zeta,
                       ~sweep(.x %*% ortho_coef$Q, 2, .y, FUN = "*") |>
                         rowSums())

  # Add mean curve to sample paths

  mu_idx <- lapply(seq_along(points_list), function(i) {
    sapply(points_list[[i]], function(x) which.min(abs(x - mu$t)))
  })

  paths_shift <- purrr::map2(paths, mu_idx, ~.x + mu$x[.y])

  sigma_list <- purrr::map(points_list,
                           ~rnorm(n = length(.x),
                                  sd = sigma_het(sigma0 = sigma0,
                                                 hetero = TRUE,
                                                 t_vec = .x)))

  # Return relevant objects
  lapply(seq_along(paths_shift), function(x) {
    list(curves = list(t = points_list[[x]],
                       x = paths_shift[[x]] + sigma_list[[x]]),
         mu = list(t = mu$t[mu_idx[[x]]],
                   x = mu$x[mu_idx[[x]]]),
         zeta = zeta[[x]],
         v = fourier_base_grid[[x]] %*% ortho_coef$Q,
         sigma = sigma_list[[x]])
  })

}

#' Computes the covariance function for some input eigenvalues
#' and eigenfunctions using Mercer's theorem.
#'
#' @param x_val Vector of eigenvalues.
#' @param y_fun Matrix, with each column containing the j-th eigenfunction
#' @return Matrix, containing the K x K covariance function, where K
#' is the number of rows of `y_fun`.
#' @export

cov_mercer <- function(x_val, y_fun) {

  val_vec <- sweep(y_fun,
                   MARGIN = 2,
                   x_val,
                   FUN = "*")

  sapply(seq_along(x_val),
         function(j) tcrossprod(val_vec[, j]), simplify = "array") |>
    apply(MARGIN = c(1, 2), sum)

}


