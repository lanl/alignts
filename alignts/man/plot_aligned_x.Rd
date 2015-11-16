\name{plot_aligned_x}
\alias{plot_aligned_x}
\title{ Plot observed time series output aligned to the latent profile}
\description{
Displays the series observations values aligned in latent time to the latent characteristic profile. Observed time series values are also scaled by the inverse of the scaling values that take the latent curve to the data. 
}
\usage{
plot_aligned_x(t, x,aligned_obj,type="p")
}
\arguments{
  \item{t}{ list of time indices for each replicate time series }
  \item{x}{ list of output values for each replicate time series }
  \item{aligned_obj}{ alignment model output from \code{align_series_EM} or \code{predict_series_EM} }
  \item{type}{ plot type. Defaults to "p" to avoid sloppy lines when time series has missing observations. "l" is preferred when the series are fully observed. }
}
\value{None.}
\author{
M.J. Grosskopf \email{mike.grosskopf@gmail.com}
}
\note{}

\seealso{
\code{\link{plot_raw}},\code{\link{plot_residuals}},\code{\link{warp_example}}
}
\examples{
## load example data
data(uvDiode)

## Set up a vector of time indices
##     Function supports a vector in place of a list if all
##     observation lists are the same length
t <- 1:nrow(uvDiode)

## Fit the time series alignment model
aligned <- align_series_EM(t, uvDiode,lambda=0.1,iters=50,upsample_factor = 1.2,tol=1.e-5)
  
## Plot aligned time series from the above model
plot_aligned_x(t,uvDiode,aligned)
}
