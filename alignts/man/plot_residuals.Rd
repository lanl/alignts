\name{plot_residuals}
\alias{plot_residuals}
\title{ Plot residuals for each time series after application of the alignment model}
\description{
Streamlines the display of the residuals of the model. 
}
\usage{
plot_residuals(t, residuals)
}
\arguments{
  \item{t}{ list of time indices for each replicate time series }
  \item{residuals}{ list of residual values for each replicate time series }
}
\value{None.}
\author{
M.J. Grosskopf \email{mike.grosskopf@gmail.com}
}
\note{}

\seealso{
\code{\link{plot_aligned_x}},\code{\link{plot_raw}}
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

## Obtain residuals from the time series alignment model
residuals <- get_residuals(t,uvDiode,aligned,standardize=TRUE)

## Plot the residuals acquired above
plot_residuals(t,residuals)

}
