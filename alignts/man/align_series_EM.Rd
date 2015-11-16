\name{align_series_EM}
\alias{align_series_EM}
\title{Align Time-warped, Scaled, Time-Series to a Characteristic Latent Profile}

\description{
  A hidden Markov Model (HMM) approach to taking replicate time series that have a characteristic behavior, but may be warped in time and multiplicatively scaled in output value. 
}

\usage{
align_series_EM(t, x, J=5, upsample_factor=2,buffer=NULL, tol=1.e-5, transition_tau=NULL, transition_scale=c(0.5,0.25), noise=NULL, init_t_range=NULL, scales=1, iters=50,lambda=0,parallel=TRUE,periodic=FALSE,get_viterbi=TRUE)
}

\arguments{
  \item{t}{ list of time indices for each replicate time series }
  \item{x}{ list of output values for each replicate time series }
  \item{J}{ number of allowed states in latent time to which the HMM can transition (default=5) }
    \item{init_t_range}{ scalar setting the number of time states accessible by the first observation of a series. Allows initial location of a series to vary without necessitating large values of \code{J} for the temporal transitions }
  \item{upsample_factor}{ multiplicative factor for setting the number of latent time states. If result is not an integer the value is floored. Recommended by Listgarten (2004) to be 2 (default=2.0) }
  \item{buffer}{ number of additional time states to be added to either end of the latent profile to fully set the length. The total number of time states of the latent series is \code{floor}(\code{upsample_factor}*\code{max}(\code{lengths}(\code{x})))+2*\code{buffer}. The default for \code{buffer} is 5\% of the longest replicate series as recommended in Listgarten (2004) }
  \item{tol}{ covergence tolerance on the penalized log-likelihood. Convergence determined when the fractional change in penalized log-likelihood is less than \code{tol}. Default is 1.e-5.}
  \item{transition_tau}{ vector of transition probabilities to the \code{J} states ahead in time from the current state. Default is to set all \code{J} transitions to have probability 1/\code{J} }
  \item{transition_scale}{ length 2 vector containing the transition probabilites on the scale states. The first value contains the probability of staying at the same scale, the second gives the probability of transitioning out to an adjacent scale state. The probability of transitioning up 1 and down 1 are both equal to the second value. (default = (0.5,0.25)) }
  \item{noise}{ scalar to initialize the noise variance for the E-M algorithm. Defaults to the square of 15\% of the range of the output values, as recommended in Listgarten (2004) }
  \item{scales}{ vector containing the values for the latent scaling parameter. Defaults to one scale state with value of 1. }
  \item{iters}{ scalar maximum number of iterations for the E-M algorithm. (default = 50) }
  \item{lambda}{ smoothness parameter on the latent trace. (default = 0.) } 
  \item{n_times}{ number of times in observed series if fully observed. Setting to NULL allows it to be set by the number of unique time values in the data. (default = NULL) } 
  \item{parallel}{ logical indicating whether OpenMP parallelization should be used where available. (default = TRUE) }
  \item{periodic}{ logical indicating whether learned trace should allow transitions from end of path to beginning. Work in progress. (default = FALSE) }
  \item{get_viterbi}{ logical indicating whether the maximum probability path for each series should be returned. This path is found using the Viterbi algorithm. (default = TRUE) }
}
\value{
  \item{ z }{ vector of latent profile from the alignment procedure }
  \item{ tau }{ vector of latent times corresponding the to the discrete sampling of the latent profile }
  \item{ u }{ vector of global scaling parameters for the alignment model }
  \item{ noise }{ scalar error variance estimate from the alignment model }
  \item{ scales }{ vector of values for the scale states of the model }
  \item{ states }{ matrix of the scale and time state index pair for each HMM state }
  \item{ viterbi_path }{ matrix where the jth column is the state indices of the maximum probability sequence for the jth observed series }
  \item{ viterbi_z }{ matrix where the jth column is the latent profile vector of the maximum probability sequence for the jth observed series }
  \item{ viterbi_tau }{ matrix where the jth column is the latent time vector of the maximum probability sequence for the jth observed series }
  \item{lambda}{ scalar penalty parameter }
  \item{log_like}{ the penalized log-likelihood at each iteration of the E-M algorithm }
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
}
\author{
M.J. Grosskopf <mike.grosskopf@gmail.com>
}

\references{
Listgarten, J., Neal, R. M., Roweis, S. T., & Emili, A. (2004). Multiple alignment of continuous time series. 
In Advances in neural information processing systems (pp. 817-824).
}
