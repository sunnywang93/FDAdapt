% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameters.R
\name{smooth_curves}
\alias{smooth_curves}
\title{Smooth curves with a given bandwidth}
\usage{
smooth_curves(data, grid, bandwidth)
}
\arguments{
\item{data}{List of curves, with each element/curve containing two entries:
\itemize{
\item \strong{$t} Vector of time points along each curve.
\item \strong{$x} Vector of observed points along each curve.
}}

\item{grid}{Vector / List containing the grid of points to smooth curves. If
the input is a list, then smoothing is performed matching the i-th curve
in \code{data} and the i-th element on the grid.}

\item{bandwidth}{Numeric containing the bandwidth.}
}
\value{
List, containing the smoothed curves.
}
\description{
Performs smoothing of curves for a given bandwidth using the
Nadaraya-Watson smoother.
}
