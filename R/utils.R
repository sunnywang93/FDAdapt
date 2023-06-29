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





