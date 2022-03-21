\name{warp_steps_illustrate}
\alias{warp_steps_illustrate}
\title{ Display the alignment model in steps, adding one alignment parameter at each step }
\description{
Display the alignment broken down by alignment parameter - the time-warping, the local scaling, and the global scaling - for describing the model and diagnosing possible model fit issues.
}
\usage{
warp_steps_illustrate(t, x, aligned_obj, series_ind, plot_every=1,plot_separate = FALSE)
}
\arguments{
  \item{t}{ list of time indices for each replicate time series }
  \item{x}{ list of output values for each replicate time series }
  \item{aligned_obj}{ alignment model output from \code{align_series_EM} or \code{predict_series_EM} }
  \item{series_ind}{ index of the replicate time series to plot }
  \item{plot_every}{ Allows subsetting of total series observations in the plot for rendering time and visual clarity. Default of 1 plots every observed value in the series. }
  \item{plot_separate}{ logical to indicate whether the plots should be tiled on one window or plotted in separate windows }
}
\value{None.}
\author{
M.J. Grosskopf \email{mike.grosskopf@gmail.com}
}
\note{}

\seealso{
\code{\link{warp_example}}
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
  
## Plot the sixth replicate series, showing the point on the latent curve that
##     each observed value aligned to, broken down by each alignment parameter
warp_steps_illustrate(t,uvDiode,aligned,6)
}
