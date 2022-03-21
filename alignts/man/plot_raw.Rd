\name{plot_raw}
\alias{plot_raw}
\title{ Plot raw replicate time series }
\description{
Streamlines creation line plots of the replicate time series for quick visualization
}
\usage{
plot_raw(t, x, lwd=2, xlab="", ylab="")
}
\arguments{
  \item{t}{ list of time indices for each replicate time series }
  \item{x}{ list of output values for each replicate time series }
  \item{lwd}{ line width of the time series (default=2) }
  \item{xlab}{ x-axis label on the plot (default="") }
  \item{ylab}{ y-axis label on the plot (default="") }
}
\value{None.}
\author{
M.J. Grosskopf \email{mike.grosskopf@gmail.com}
}
\note{}

\seealso{
\code{\link{plot_aligned_x}},\code{\link{plot_residuals}}
}
\examples{
## load example data
data(uvDiode)

## Set up a vector of time indices
##     Function supports a vector in place of a list if all
##     observation lists are the same length
t <- 1:nrow(uvDiode)

## Plot replicate time-series
plot_raw(t,uvDiode)
}
