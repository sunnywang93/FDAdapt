% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FPCA.R
\name{FPCA}
\alias{FPCA}
\title{Performs functional principal components analysis}
\usage{
FPCA(
  data,
  grid_smooth,
  grid_bw,
  grid_param,
  cv_set,
  quantile,
  inflate_bw = TRUE,
  center = TRUE,
  interp_method = "linear",
  nelements = 10,
  gamma_H,
  gamma_L,
  n_knots,
  h_power,
  intp_param,
  true_params,
  diag_cor = TRUE
)
}
\arguments{
\item{data}{List of curves, with each element/curve containing two entries:
\itemize{
\item \strong{$t} Vector of time points along each curve.
\item \strong{$x} Vector of observed points along each curve.
}}

\item{grid_smooth}{Vector of sampling points to smooth curves. Eigenfunctions
will be computed at these sampling points.}

\item{grid_bw}{Vector of points for bandwidth optimization.}

\item{grid_param}{Vector of sampling points to estimate parameters.}

\item{cv_set}{Numeric, containing the number of curves used for learning
the cross-validation bandwidth in presmoothing.}

\item{quantile}{Numeric, containing the quantile selected from the set of
cross-validation bandwidths.}

\item{inflate_bw}{Boolean, indicating whether to correct the bandwidth due
to discretization error from estimating the regularity.}

\item{center}{Boolean, indicating whether to center curves.}

\item{nelements}{Numeric, indicating the number of eigen-elements to keep.}

\item{gamma_H}{Numeric, indicating the gamma to be used for estimation of
\code{H}. See \code{estimate_regularity}.}

\item{gamma_L}{Numeric, indicating the gamma to be used for estimation of
\code{L}. See \code{estimate_regularity}.}

\item{n_knots}{Numeric, number of knots used for smoothing H and L, where
splines are used.}

\item{h_power}{Numeric, power to raise the bandwidth to when estimating
sigma for the purposes of diagonal correction.}

\item{intp_param}{Boolean, where \code{TRUE} indicates using presmoothing + interpolation
to estimate regularity parameters, as compared to only presmoothing.}

\item{interp_type}{String, indicating the type of interpolation to perform
for the parameters from \code{grid_param} onto \code{grid_smooth}. Options include
c("linear", "constant", "nearest", "spline", "cubic").}
}
\value{
List, containing the following elements:
\itemize{
\item \strong{$params} List containing the estimated parameters.
\item \strong{$bw} Numeric containing the bandwidth used for smoothing curves.
\item \strong{$evalues} Vector containing the normalised eigenvalues.
\item \strong{$efunctions} Matrix containing the normalised eigenfunctions, with the j-th
column representing the j-th eigenfunction.
}
}
\description{
Perform FPCA of functional data adaptively using the local regularity of
curves. See references for more details on the algorith.
}
