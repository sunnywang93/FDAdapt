#' Generate exponential grid for optimisation
#'
#' @param from Starting point of grid.
#' @param to End point of grid.
#' @param length.out Number of points in the grid.
#' @returns A vector containing grid values.
#' @examples
#' lseq(from = 1, to = 100, length.out = 51)
#' @export

lseq <- function (from = 1, to = 100, length.out = 51) {
  exp(seq(log(from), log(to), length.out = length.out))
}

#' @export
epa_kernel <- function(y) {
  3/4 * (1 - y^2) * (abs(y) <= 1)
}


#' @export
epa_shrinkage <- function(u){
  ifelse(abs(u) <= 1, 0.75*(1 - u**2), .Machine$double.eps)
}



#' @export
bertin_kernel <- function(x, beta){
  ((1 + beta)/ (2 * beta))  * (1 - abs(x)^beta) * (abs(x) <= 1)
}



ISE <- function(t, x, y, trimmed = FALSE){
  if (trimmed) {
    cond <- (t > 0.05) & (t < 0.95)
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

#' @export
list2cai <- function(data){
  seq_along(data) |>
    lapply(function(idx) {
      data.frame(obs = idx, time = data[[idx]]$t, x = data[[idx]]$x)
    }) |>
    (\(x) do.call("rbind", x))()
}

#' @export
mean_ss <- function(curves, grid){
  curves_ <- list2cai(curves)
  mod <- stats::smooth.spline(curves_$time, curves_$x)
  stats::predict(mod, grid)$y
}

#' @export
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


# generate_mean_curve <- function(k_length, grid_t = seq(0, 1, length.out = 101),
#                                 alpha, shift = 0, scale_mu = 1) {
#   Z <- runif(k_length, -sqrt(3), sqrt(3))
#   #Z <- rnorm(k_length)
#   xi <- c()
#   # for(i in 1:k_length) {
#   #   xi_k <- (-1)^(i+1)*i^(-alpha)
#   #   xi <- c(xi, xi_k)
#   # }
#   for(i in 1:k_length) {
#     xi_k <- 1/((i - 1/2)^alpha*pi^2)
#     xi <- c(xi, xi_k)
#   }
#   k <- seq(1, k_length, by = 1)
#
#   # mu_kt <- sapply(grid_t, function(t) {
#   #   Z * sqrt(xi) * cos(k * pi * t)
#   # })
#
#   mu_kt <- sapply(grid_t, function(t) {
#     sqrt(2) * Z * sqrt(xi) * sin((k - 0.5) * pi * t)
#   })
#
#
#   mu_t <- scale_mu * colSums(mu_kt) + rep(shift, length(grid))
#   tibble(t = grid_t, mu = mu_t)
# }


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


#' @export
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

#' @export
normalise_eigen <- function(covariance, nelements = 10) {
  eelements <- eigen(covariance, symmetric = TRUE)
  evalues <- eelements$values[seq(nelements)]
  efunctions <- eelements$vectors[, seq(nelements)]
  evalues_norm <- evalues / nrow(covariance)
  efunctions_norm <- sapply(seq(nelements), function(j) efunctions[, j] *
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


#' Perform FPCA by interpolation
#'
#' Linearly interpolate curves and obtain eigen-elements by eigen-decomposition
#' of the empirical mean and covariance functions
#'
#' @param curves List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param grid_smooth Vector of sampling points to perform interpolation.
#' @param nelements Numeric containing the number of eigen-elements to return.
#' @returns List, with the following components
#' - **evalues** Vector, containing the normalised eigenvalues.
#' - **efunctions** Matrix, containing the normalised eigenfunctions.
#' @export
eigen_interpolate <- function(curves, grid_smooth, nelements) {

  smoothed_curves <- sapply(curves, function(i) {
    stats::approx(x = i$t,
                  y = i$x,
                  xout = grid_smooth,
                  method = "linear",
                  yleft = i$x[1],
                  yright = i$x[length(i$x)])$y
  })

  mean <- rowMeans(smoothed_curves, na.rm = TRUE)
  centered_curves <- smoothed_curves - mean
  cov_all <- lapply(seq_along(curves),
                    function(i) tcrossprod(centered_curves[, i]))

  cov <- Reduce('+', cov_all) / length(cov_all)

  eelements <- normalise_eigen(cov, nelements)
  list(evalues = eelements$values, efunctions = eelements$vectors)

}

#' Normalise the sign of estimated eigenfunctions
#'
#' Since estimated eigenfunctions are only defined up to a sign, normalisation
#' is required to have identical signs with the true eigenfunctions for
#' the purposes of comparison.
#'
#' @param efunction Matrix containing the estimated eigenfunctions.
#' @param efunction_true Matrix containing the true eigenfunctions.
#' @return Matrix, containing the normalised eigenfunctions.
#' @export
normalise_sign <- function(efunction, efunction_true) {
  sapply(seq(ncol(efunction)), function(j) {
    c(sign(t(efunction[, j]) %*% efunction_true[, j])) * efunction[, j]
  })
}


#' Perform one dimensional interpolation
#'
#' Interpolation using `approx` from `stats`, where the extrapolation is done
#' by using the closest observed point.
#'
#' @param x Vector, providing the sampling points to be interpolated from.
#' @param y Vector, providing the observed points to be interpolated from.
#' @param xout Vector, providing where interpolation should take place.
#' @param type String, specifying the interpolation method. Choices are either
#' "linear" or "constant".
#' @returns List, with the following components
#' - **x** Vector, containing the interpolated sampling points.
#' - **y** Vector, containing the interpolated points.
#' @references Wang S., Patilea V., Klutchnikoff N. (2023+) - Adaptive Functional
#' Principal Components Analysis
#' @export

interpolate1D <- function(x, y, xout, type = "linear") {

  stats::approx(x = x,
                y = y,
                xout = xout,
                method = type,
                yleft = y[1],
                yright = y[length(y)],
                ties = mean)

}


#' Perform two dimensional interpolation
#'
#' Bilinear interpolation using `intp` package.
#'
#' @param x Vector, providing the `x` coordinates to be interpolated from.
#' @param y Vector, providing the `y` coordinates to be interpolated from.
#' @param z Matrix, containing the `z(x, y)` values to be interpolated from.
#' @param xout Vector, providing the `x` coordinates where interpolation should
#' take place.
#' @param yout Vector, providing the `y` coordinates where interpolation should
#' take place.
#' @returns Matrix, containing the interpolated points for every `z(xout, yout)`.
#' @references Wang S., Patilea V., Klutchnikoff N. (2023+) - Adaptive Functional
#' Principal Components Analysis
#' @export

interpolate2D <- function(x, y, z, xout, yout) {

  zout <- interp::bilinear.grid(x = y,
                                y = y,
                                z = z,
                                nx = length(xout),
                                ny = length(yout))

  zout$z
}


#' Estimates smoothed mean from a dataset
#'
#' Smooths the empirical mean of a dataset using a lasso penalty.
#' @param df Dataframe, with rows indexing the curves and columns indexing
#' the time points.
#' @param k Number of basis functions to be used in fit
#' @returns A `glmnet` object, which can be used as input to `predict_mean` to
#' evaluate the smoothed mean on a specified grid.
#' @export

learn_mean <- function (df, k = 50) {
  true_mu <- unname(colMeans(df, na.rm = TRUE))
  m <- ncol(df)
  t <- seq(0, 1, length.out = m)
  tfeatures <- matrix(NA, ncol = 2 * k + 1, nrow = length(t))
  tfeatures[, 1] <- t
  for (j in 1:k) {
    tfeatures[, j + 1] <- sqrt(2) * cos(2 * j * pi * t)
    tfeatures[, 2 * k + 2 - j] <- sqrt(2) * sin(2 * j * pi * t)
  }
  glmnet::glmnet(x = tfeatures, y = true_mu, alpha = 1)
}


#' Evaluates the smoothed mean on a grid
#'
#' Given a `glmnet` object, for example from the `learn_mean` function,
#' `predict_mean` computes the function at specific evaluation points.
#' @param u Vector of sampling points on which the mean which be computed at.
#' @param model `glmnet` object containing the mean curve.
#' @param lambda Numeric, representing the value of the regularisation parameter.
#' @param k Number of basis functions used to fit the data.
#' @param scale Boolean. If `TRUE`, scales the output by the global mean.
#' @returns A numeric vector containing the mean curve evaluated on the vector of
#' grid points `u`.
#' @export
predict_mean <- function (u, model, lambda, k = 50, scale) {
  features <- matrix(NA, ncol = 2 * k + 1, nrow = length(u))
  features[, 1] <- u
  for (j in 1:k) {
    features[, j + 1] <- sqrt(2) * cos(2 * j * pi * u)
    features[, 2 * k + 2 - j] <- sqrt(2) * sin(2 * j * pi *
                                                 u)
  }
  if(scale) {
    muhat <- glmnet::predict.glmnet(model, newx = features, s = lambda)[, 1]
    muhat - mean(muhat)
  } else {
    glmnet::predict.glmnet(model, newx = features, s = lambda)[, 1]
  }
}

#' Check if curves are sampled on a common design scheme
#'
#' @export
check_common <- function(x, y) {

  if(identical(x, y)) {
    y
  } else {
    FALSE
  }

}


#' Computes the maximum (over observed points) Nadaraya-Watson weights
#'
#' @param curves List, containing the following elements:
#' -**t** Vector, containing the sampling points of curves.
#' -**x** Vector, containing the observed points.
#' @param x Vector, containing the evaluation points.
#' @param bw Vector, containing the grid of bandwidth points.
#' @returns Array of dimension `G x N x H`, corresponding to the
#' number of evaluation points, curves, and bandwidth points respectively.
#' @export
NW_max <- function(curves, x, bw) {

  weights_max <- lapply(curves, function(i) {
    epa_kernel(outer(outer(i$t, x, FUN = "-"),
                     bw, FUN = "/")) |>
      apply(c(2, 3), max)
  })

  weights_denom <- lapply(curves, function(i) {
    epa_kernel(outer(outer(i$t, x, FUN = "-"),
                     bw, FUN = "/")) |>
      colSums()
  })

 purrr::map2(weights_max, weights_denom, ~(.x / .y)) |>
   rapply(f = function(x) ifelse(is.nan(x), 0, x), how = "replace") |>
   abind::abind(along = 3) |>
   aperm(c(1, 3, 2))

}


#' Computes the weights for curve selection
#'
#' Given a list of functional data and a grid of points `x`, indicator functions
#' are computed which checks if there is at least one point along a curve
#' in the neighbourhood (defined by a bandwidth) of an evaluation point in `x`.
#'
#' @param curves List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param x Vector, containing the evaluation points.
#' @param h Vector or numeric, containing the bandwidth.
#' @returns Array of dimension `G x N x H`, corresponding to the
#' number of evaluation points, curves, and bandwidth points respectively.
#' @export
#'
curves_select <- function(curves, x, h) {

  wi <- lapply(curves, function(curve) {
    abs(outer(curve$t, x, FUN = "-")) |>
      outer(h, FUN = "<=") |>
      (\(x) (colSums(x) >= 1) * 1)()
    }) |>
    abind::abind(along = 3) |>
    aperm(c(1, 3, 2))

  if(dim(wi)[3] == 1) {
    wi[,, 1]
  } else {
    wi
  }

}

#' Computes part of bias term involving kernels in adaptive estimation
#'
#' @param H Vector containing the Hölder exponents.
#' @param L Vector containing the Hölder constants.
#' @param cst Vector containing the constant term depending on the kernel.
#' @param h_grid Vector containing the grid of bandwidths.
#' @returns Vector with the same length as `H`, `L` and `cst`.
#' @export

bias_kernel <- function(H, L, cst, h_grid) {

  identicalValue <- function(x,y) if (identical(x,y)) x else FALSE

  if(is.logical(Reduce(identicalValue, length(H), length(L), length(cst)))) {
    stop("H, L and cst must be of the same length!")
  }

  outer(h_grid, 2 * H, FUN = "^") |>
    sweep(MARGIN = 2, STATS = L**2, FUN = "*") |>
    sweep(MARGIN = 2, STATS = cst, FUN = "*")

}

#' Interpolates elements of a list
#'
#' Performs 1D and 2D interpolation of elements in a list
#'
#' @param list_in List, containing
#' -**t** Vector of evaluation points in which the elements of the list were
#' computed on.
#' -Vector or matrices of to be interpolated.
#' @param xout Vector, containing the desired output evaluation points.
#' @returns List, with each element interpolated onto the grid of `xout`.
#' @export

intp_list <- function(list_in, xout) {

  # Interpolate 1D-parameters onto smoothing grid
  if(length(list_in$t) != length(xout)) {

    list_sub <- list_in[purrr::map_lgl(list_in, ~is.vector(.x) && length(.x) > 1)]

    list_smooth <- lapply(within(list_sub, rm(t)),
                          function(i) interpolate1D(x = list_in$t,
                                                    y = i,
                                                    xout = xout)$y)

    # Interpolate 2D-parameters
    list2D_smooth <- lapply(list_in[purrr::map_lgl(list_in, ~is.matrix(.x))],
                            function(X) interpolate2D(x = list_in$t,
                                                      y = list_in$t,
                                                      z = X,
                                                      xout = xout,
                                                      yout = xout)
    )


    # Append to list of parameters
    c(list(t = xout),
      list_smooth,
      list_in[purrr::map_lgl(list_in, ~length(.x) == 1)],
      list2D_smooth
      )

  } else {
    list_in
  }


}

#' Computes the variance rate term of the covariance function
#'
#' This quantity computes the term that governs that rate of convergence of
#' the variance term for the covariance function.
#'
#' @param curves List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param x Vector containing the evaluation points in one dimension. The final
#' quantity is computed over the cartesian product of (`x`, `x`).
#' @param bw Vector containing the grid of bandwidths.
#' @returns List, containing two elements:
#' - **WN_bi** Array of dimension G x G x H, where G is the length of `x` and
#' `H` is the length of `bw`. This term contains the effective sample size for
#' a given bandwidth in the estimation of the covariance function.
#' - **Ngamma** Array of dimension G x G x H, where G is the length of `x` and
#' `H` is the length of `bw`. This term contains the variance term in the risk
#' bounds of the covariance function, which contributes to the rate of convergence.
#' @export

N_gamma <- function(curves, x, bw) {

  wi <- curves_select(curves = curves,
                      x = x,
                      h = bw)

  WN_bi <- sapply(seq_along(curves), function(i) {
    sapply(seq_along(bw), function(h) {
      tcrossprod(wi[,i,h])
    }, simplify = "array")
  }, simplify = "array") |>
    aperm(c(4, 1, 2, 3)) |>
    colSums()

  Wm_max <- NW_max(curves = curves,
                   x = x,
                   bw = bw)

  Ngamma <- sapply(seq_along(curves), function(i) {
    sapply(seq_along(bw), function(h) {
      tcrossprod(wi[,i,h] * Wm_max[,i,h], wi[,i,h])
    }, simplify = "array")
  }, simplify = "array") |>
    aperm(c(4, 1, 2, 3)) |>
    colSums() |>
    (\(x) (x / WN_bi**2)**(-1))()

  Ngamma[is.nan(Ngamma)] <- 0

  list(
    WN_bi = WN_bi,
    Ngamma = Ngamma
  )

}


#' Computes the diagonal bias of the covariance function
#'
#' @param curves List of curves, with each element/curve containing two entries:
#' - **$t** Vector of time points along each curve.
#' - **$x** Vector of observed points along each curve.
#' @param t_grid Vector of evaluation points where curves are smoothed.
#' @param h Numeric containing the bandwidth (usually an optimal one).
#' @param zeta Numeric containing the power of the bandwidth.
#' @param wm_list List, where each element contains an Mi x G matrix of
#' Nadaraya-Watson weights, where Mi is the number of points along curve i and
#' G is the length of `t_grid`.
#' @param W Matrix of dimension G x N, where N is the number of curves. It
#' contains the weights for curves selection, resulting from the function `curves_select`.
#' @param WN Matrix of dimension G x G, containing the number of effective sample
#' size used for covariance function estimation.
#' @returns Matrix, containing the diagonal bias.
#' @export

diagonal_bias <- function(curves, t_grid, h, zeta, wm_list, W, WN) {


  sigma_diag <- estimate_sigma(data = curves,
                               sigma_grid = t_grid,
                               h = h,
                               h_power = zeta)


  diag_sum <- purrr::imap(wm_list,
                          ~crossprod(.x) * tcrossprod(W[, .y])) |>
    (\(x) Reduce('+', x))()

  tcrossprod(sigma_diag) * diag_sum / WN

}




#' Computes the empirical mean of functional data
#'
#' Given a learning and online set of curves, the empirical mean on the
#' sampling points of the online set is computed. The empirical mean is first
#' computed on the commonly sampled points of the learning set (assumed to be
#' densely observed, without noise), before interpolation is performed onto
#' the sampling points of the online set (assumed to be on a randomly designed
#' grid).
#'
#'
#' @param curves_learn List of learning curves, containing the following as elements:
#' - **$t** Vector of sampling points.
#' - **$x** Vector of observed points.
#' @param curves_online List of online curves, containing the following as elements:
#' - **t** Vector of sampling points.
#' - **$x** Vector of observed points.
#' @returns List, containing the mean function at the sampling points of the
#' online set and the interpolated points from the learning set.
#' @export

mean_emp <- function(curves_learn, curves_online) {

  if(is.logical(Reduce(check_common, purrr::map(curves_learn, ~.x$t)))) {
    stop("learning curves have must common sampling points!")
  }

  mu_learn <- Reduce('+', purrr::map(curves_learn, ~.x$x)) / length(curves)

  mu_online <- purrr::map(curves_online,
                          ~interpolate1D(x = curves_learn[[1]]$t,
                                         y = mu_learn,
                                         xout = .x$t))
  list(
    mu_learn = list(t = curves[[1]]$t,
                    x = mu_learn),
    mu_online = purrr::map(mu_online, ~list(t = .x$x,
                                            x = .x$y))

  )

}

#' Computes the empirical covariance function for functional data
#'
#' Given a list of commonly sampled functional data, the empirical mean
#' is computed on the sampling points, before interpolation if an optional grid
#' of output is specified. Note that if the evaluation points are specified,
#' the covariance function is first constructed by centering each curve at
#' the observed points, before interpolation is done.
#'
#' @param curves_learn List of learning curves, containing the following as elements:
#' - **$t** Vector of sampling points.
#' - **$x** Vector of observed points.
#' @param curves_online List of online curves, containing the following as elements:
#' - **t** Vector of sampling points.
#' - **$x** Vector of observed points.
#' @returns List, containing the following elements:
#' -**$cov_learn** Matrix containing the covariance function from the learning set.
#' -**$cov_online** List, containing the sampling points from the online set and
#' the interpolated covariance function.
#' @export

cov_emp <- function(curves_learn, curves_online) {

  if(is.logical(Reduce(check_common, purrr::map(curves_learn, ~.x$t)))) {
    stop("learning curves have must common sampling points!")
  }

  mu_emp <- mean_emp(curves_learn = curves_learn,
                     curves_online = curves_online)

  cov_learn <- purrr::map(curves_learn, ~.x$x - mu_emp$mu_learn$x) |>
    purrr::map(~tcrossprod(.x)) |>
    (\(x) Reduce('+', x) / length(x))()

  cov_out <- purrr::map(curves_online,
                        ~interpolate2D(x = curves_learn[[1]]$t,
                                       y = curves_learn[[1]]$t,
                                       z = cov_learn,
                                       xout = .x$t,
                                       yout = .x$t)
                        )

  list(cov_learn = cov_learn,
       cov_online = purrr::map2(curves_online, cov_out,
                                  ~list(t = .x$t, cov = .y)),
       mu_learn = mu_emp$mu_learn,
       mu_online = mu_emp$mu_online
       )


}


#' Computes PACE estimates based on a learning and online set of curves
#'
#' Given a list of commonly sampled functional data, the conditional
#' expectation of the scores given the observed values are computed. The
#' plug-in estimates of the conditional expectation are the empirical
#' mean, covariance, and eigen-elements.
#'
#' @param curves_learn List of learning curves, containing the following as elements:
#' - **$t** Vector of sampling points.
#' - **$x** Vector of observed points.
#' @param curves_online List of online curves, containing the following as elements:
#' - **t** Vector of sampling points.
#' - **$x** Vector of observed points.
#' @param sigma Numeric, containing the noise of the online curves.
#' @returns List, containing the scores for each online curve.
#' @export

PACE_emp <- function(curves_learn, curves_online, sigma) {

  if(is.logical(Reduce(check_common, purrr::map(curves_learn, ~.x$t)))) {
    stop("learning curves have must common sampling points!")
  }

  cov <- cov_emp(curves_learn = curves_learn,
                 curves_online = curves_online)

  cov_noisy <- purrr::map(cov$cov_online,
                          ~.x$cov + diag(sigma**2,
                                         nrow = nrow(.x$cov),
                                         ncol = ncol(.x$cov))
                          )

  eelements <- eigen(cov$cov_learn, symmetric = TRUE)

  efunctions_out <- purrr::map(curves_online,
                               ~apply(eelements$vectors, 2,
                                      function(v) {
                                        interpolate1D(x = curves_learn[[1]]$t,
                                                      y = v,
                                                      xout = .x$t)$y
                                        }
                                      )
                               )


  curves_cen <- purrr::map2(curves_online, cov$mu_online,
                            ~.x$x - .y$x)

  lapply(seq_along(curves_online), function(id) {
    sweep(efunctions_out[[id]], 2, eelements$values, FUN = "*") |>
      apply(2, function(v) t(v) %*% solve(cov_noisy[[id]]) %*% curves_cen[[id]])
  })


}


#' Orthonormalise the eigenfunction matrix
#'
#' Given a matrix where the columns are composed of the linearly
#' independent eigenfunctions, the QR decomposition is returned. Normalisation
#' is taking performed to make the eigenfunctions of unit norm, and have the
#' same sign as the input eigenfunctions.
#'
#' @param X Matrix, containing the eigenfunctions in the columns to be
#' orthonormalised.
#' @param t Vector containing the evaluation points of the eigenfunction.
#' @return Matrix, containing the orthonormalised eigenfunctions.
#' @export

ortho_funmat <- function(X, t) {

  if(nrow(X) != length(t)) {
    stop("Number of evaluation points do not match in t and X!")
  }

  Q <- qr.Q(qr(X))

  norm_cst <- apply(Q, 2, function(x) sqrt(pracma::trapz(t, x**2)))

  X_ortho <- sweep(Q, 2, norm_cst, FUN = "/")

  normalise_sign(efunction = X_ortho,
                 efunction_true = X)

}




#' Powerconsumption data set
#'
#' A subset of the powerconsumption dataset from the UCI Repository
#' @format A tibble with 1358 rows and 1440 columns, with rows representing
#' the daily curves and columns representing 1 minute samples over a 24-hour
#' period.
#' \describe{
#' Measurements of voltage consumption in one household with a one-minute
#' sampling rate in a 24 hour period, with data taken over almost 4 years. Days with missing
#' values were removed from the dataset.
#' }
#' @source \url{https://archive.ics.uci.edu/ml/datasets/Individual+household+electric+power+consumption}
"powerconsumption"





