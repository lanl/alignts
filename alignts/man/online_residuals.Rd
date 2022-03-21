\name{online_residuals}
\alias{online_residuals}
\title{Obtaining residuals alignment of a new replicate series }

\description{
Obtains the residuals for a new series, given the alignment results from \code{online_align}. 
}

\usage{
online_residuals(t, x, aligned_obj, standardize = FALSE)
}

\arguments{
  \item{t}{ list of time indices for each replicate time series }
  \item{x}{ list of output values for each replicate time series }
  \item{aligned_obj}{  alignment model output from \code{align_series_EM} or \code{predict_series_EM} }
  \item{standardize}{ logical to determine if the residuals should be standardized by the estimated model standard error }
}
\value{
List of vectors containing the residuals for each time series.
}
\seealso{
\code{\link{online_align}}, \code{\link{get_residuals}}
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

## Obtain standardized residuals from the alignment model
residuals <- online_residuals(t,uvDiode[,1],aligned,standardize=TRUE)
}
\author{
M.J. Grosskopf <mike.grosskopf@gmail.com>
}
