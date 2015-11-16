\name{warp_example}
\alias{warp_example}
\title{ Display the alignment of each observation in one replicate series to the latent profile. }
\description{
Visualization of the alignment of a given observed series to the latent profile.
}
\usage{
warp_example(t, x, aligned_obj, series_ind, plot_every=1)
}
\arguments{
  \item{t}{ list of time indices for each replicate time series }
  \item{x}{ list of output values for each replicate time series }
  \item{aligned_obj}{ alignment model output from \code{align_series_EM} or \code{predict_series_EM} }
  \item{series_ind}{ index of the replicate time series to plot }
  \item{plot_every}{ Allows subsetting of total series observations in the plot for rendering time and visual clarity. Default of 1 plots every observed value in the series.}
}
\value{None.}
\author{
M.J. Grosskopf \email{mike.grosskopf@gmail.com}
}
\note{}

\seealso{
\code{\link{plot_aligned_x}},\code{\link{warp_steps_illustrate}}
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
  
## Plot the sixth replicate series, showing the point on the latent curve
##     that each observed value aligned to.
warp_example(t,uvDiode,aligned,6)
}
