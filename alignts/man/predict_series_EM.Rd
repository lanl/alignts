\name{predict_series_EM}
\alias{predict_series_EM}
\title{Align New Time-Series a Previously Fit Characteristic Latent Profile}

\description{
  Running a restricted E-M algorithm to align a new series to the latent profile from a hidden Markov Model (HMM) alignment to a set of training time series. Useful for fitting model parameters with cross-validation or for obtaining residuals from newly acquired time series. 
}

\usage{
predict_series_EM(t, x, aligned_obj, tol=1.e-5, J=5, init_t_range=NULL, transition_tau=NULL, transition_scale=c(0.5,0.25), iters=50,parallel=TRUE,periodic=FALSE,get_viterbi=TRUE)
}

\arguments{
  \item{t}{ list of time indices for each replicate time series }
  \item{x}{ list of output values for each replicate time series }
  \item{aligned_obj}{ previously fit alignment model using \code{align_series_EM} }
  \item{tol}{ covergence tolerance on the penalized log-likelihood. Convergence determined when the fractional change in penalized log-likelihood is less than \code{tol}. Default is 1.e-5.}
  \item{J}{ number of allowed states in latent time to which the HMM can transition (default=5)}
  \item{init_t_range}{ scalar setting the number of time states accessible by the first observation of a series. Allows initial location of a series to vary without necessitating large values of \code{J} for the temporal transitions }
  \item{transition_tau}{ vector of transition probabilities to the \code{J} states ahead in time from the current state. Default is to set all \code{J} transitions to have probability 1/\code{J}}
  \item{transition_scale}{ length 2 vector containing the transition probabilites on the scale states. The first value contains the probability of staying at the same scale, the second gives the probability of transitioning out to an adjacent scale state. The probability of transitioning up 1 and down 1 are both equal to the second value. (default = (0.5,0.25))}
  \item{iters}{ scalar maximum number of iterations for the restricted E-M algorithm. (default = 50)}
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
}
\seealso{
\code{\link{align_series_EM}}
}
\examples{
## load example data
data(uvDiode)

## Set up a vector of time indices
##     Function supports a vector in place of a list if all
##     observation lists are the same length
t <- 1:nrow(uvDiode)

## Fit the time series alignment model
aligned <- align_series_EM(t, uvDiode[,-(1:4)],lambda=0.1,iters=50,upsample_factor = 1.2,tol=1.e-5)

## Predict the states for a new set of values, without refitting the
##     shared model parameters.
new_aligned <- predict_series_EM(t,uvDiode[,1:4],aligned,iters = 50,tol=1.e-5)
}
\author{
M.J. Grosskopf <mike.grosskopf@gmail.com>
}

\references{
Listgarten, J., Neal, R. M., Roweis, S. T., & Emili, A. (2004). Multiple alignment of continuous time series. 
In Advances in neural information processing systems (pp. 817-824).
}
