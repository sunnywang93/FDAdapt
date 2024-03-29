% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dgp.R
\name{generate_mean_curve}
\alias{generate_mean_curve}
\title{Generate a mean curve.}
\usage{
generate_mean_curve(grid = seq(0, 1, length.out = 101), K = 60, alpha = 1)
}
\arguments{
\item{grid}{Vector of sampling points to generate the mean curve.}

\item{K}{Number of basis functions to use.}

\item{alpha}{Parameter controlling the regularity of the mean curve. A higher
\code{alpha} will produce smoother mean curves.}
}
\value{
A list, containing the sampling points and observed points.
}
\description{
Generates a mean curve using a combination of sin basis functions.
}
