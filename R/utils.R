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

eigen_interpolate <- function(curves, grid_smooth, nvalues) {
  smoothed_curves <- sapply(curves, function(i) {
    pracma::interp1(x = c(0, i$t, 1),
                    y = c(i$x[1], i$x, i$x[length(i$x)]),
                    xi = grid_smooth, method = "linear")
  })
  mean <- rowMeans(smoothed_curves, na.rm = TRUE)
  centered_curves <- smoothed_curves - mean
  cov <- sapply(seq_along(curves), function(i) {
    centered_curves[, i] %*% t(centered_curves[, i])
  }, simplify = "array") |>
    apply(MARGIN = c(1, 2), mean, na.rm = TRUE)
  eelements <- eigen(cov, symmetric = TRUE)
  evalues <- eelements$values[1:nvalues] / length(grid_smooth)
  efunctions <- eelements$vectors[, 1:nvalues] * sqrt(grid_smooth)
  list(evalues = evalues, efunctions = efunctions)
}

#need to be careful with the diagonal terms
eigen_kneip <- function(curves, grid_smooth, nvalues) {
  smoothed_curves <- sapply(curves, function(i) {
    pracma::interp1(x = c(-i$t[1], i$t, 2 - i$t[length(i$t)]),
                    y = c(i$x[1], i$x, i$x[length(i$x)]),
                    xi = grid_smooth, method = "nearest")
  })
  mean <- rowMeans(smoothed_curves)
  centered_curves <- smoothed_curves - mean
  cov <- sapply(seq_along(curves), function(i) {
    centered_curves[, i] %*% t(centered_curves[, i])
  }, simplify = "array") |>
    apply(MARGIN = c(1, 2), mean, na.rm = TRUE)
  eelements <- eigen(cov, symmetric = TRUE)
  evalues <- eelements$values[1:nvalues] / length(grid_smooth)
  efunctions <- eelements$vectors[, 1:nvalues] * sqrt(grid_smooth)
  list(evalues = evalues, efunctions = efunctions)
}

#' @export
normalise_sign <- function(efunction, efunction_true) {
  sapply(seq(ncol(efunction)), function(j) {
    c(sign(t(efunction[, j]) %*% efunction_true[, j])) * efunction[, j]
  })
}


#' Corrects the diagonal band of covariance function
#'
#'
#' @param mat_input Matrix, representing the covariance function.
#' @param h Bandwidth used to smooth curve.
#' @returns Matrix, with corrected diagonal band.
#' @references Patilea V., Wang S. (2022+) - Adaptive Functional
#' Principal Components Analysis
#' @export

correct_diag <- function(mat_input, h) {

  K <- floor(1 * h * ncol(mat_input))
  mat <- (lower.tri(mat_input, diag = TRUE) * 1) * mat_input
  for(g in 2:ncol(mat)) {
    current <- mat[max(K + 1, g - 1) , 1]
    for(i in (g-1):floor((g + 1) / 2)) {
      j <- g - i
      if(abs(i - j) > K) {
        current <- mat[i, j]
      } else {
        mat[i, j] <- current
      }
    }
  }

  for(g in (2 * ncol(mat)):(ncol(mat) + 2)) {
    current <- mat[ncol(mat) , min(ncol(mat) - K, g - ncol(mat))]
    for(i in ncol(mat):floor((g + 1)/ 2)) {
      j <- g - i
      if(abs(i - j) > K) {
        current <- mat[i, j]
      } else {
        mat[i, j] <- current
      }
    }
  }
  mat + t(mat) - diag(diag(mat))
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





