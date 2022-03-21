\name{uvDiode}
\alias{uvDiode}
\docType{data}
\title{
UV diode spectrum data
}
\description{
UV diode spectrum from Listgarten (2007) paper. 21 series with 134 values in each for alignment to a common continuous series. Used for multiclass alignment by Listgarten, but for a single class example here.
}
\usage{data("uvDiode")}
\format{
  A data frame with 21 time series each with 134 observations
  \describe{
    \item{\code{series1}, \code{series2},..., \code{series21}}{a 134 row by 21 column data frame}
  }
}
\details{
}
\source{
\url{http://www.cs.toronto.edu/~jenn/CPM/}
}
\references{
Bayesian Detection of Infrequent Differences in Sets of Time Series with
Shared Structure.  Jennifer Listgarten, Radford M. Neal, Sam T. Roweis, Rachel
Puckrin and Sean Cutler, Advances in Neural Information Processing
Systems 19, MIT Press, Cambridge, MA, 2007 (NIPS 2006).
}
\examples{
data(uvDiode)

## Set up a vector of time indices
##     Function supports a vector in place of a list if all
##     observation lists are the same length
t <- 1:nrow(uvDiode)

## Plot replicate time-series
plot_raw(t,uvDiode)
}
\keyword{datasets}
