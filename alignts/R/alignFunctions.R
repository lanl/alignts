##########
# Function definitions for aligning replicate time series
#   and inferring the latent series profile 
#   based on Listgarten, Neal, Roweis, Emili
#
# Required libraries:
#   NONE
#
##########

##### 
# Function for getting the length of each member of a list
#     Exists already in R 3.2 but not before
#
lengths <- function(x){
  c(mapply(length,x))
}

#####
# Function to check to ensure that all scale and tau states
#     in x are equal to those in y. No idea why "identical"
#     was failing with this, but whatever
#
check_states <- function(x,y) all(x==y)

#####
# Randomly sample values from neighbors to upsample
#   a new series from a shorter one. For initializing
#   the latent profile. Upsample x to be of length M
#
up_sample <- function (x, M) {
  n <- length(x)
  if(M%%n==0){
    fac <- M/n
    tmp <- as.vector(t(matrix(rep(x,fac),ncol=fac)))
  }else{
    inds  <- rep(1:n,M/n+1)[1:M]
    tmp   <- x[sort(inds)]
  }
  tmp
}

#####
# Calculated the probability-weighted SSE given a state and set
#     of state probabilities
#
state_sse <- function(q, t, state_prob, x, z, tau, u, scale){
  not_na  <- which(!is.na(x))
  sum(state_prob[not_na]*(x[not_na]-z[tau[t]]*u*scale[q])^2)
}

#####
# Update sigma^2 using equation on Page 4 of paper
#
m_step_noise <- function(K,M,Q,state_prob,x,z,tau,u,scale,n_k){
  states <- expand.grid(t=1:M,q=1:Q)
  noise  <- 0
  # Make sure this is the correct N when there are NA values
  #     in the observations. It could be that this should be
  #     different.
  N <- sum(n_k)
  for (k in 1:K){
    noise <- noise + sum(sapply(1:nrow(states),
                                function(i){state_sse(states$q[i],states$t[i], state_prob[[k]][i,], 
                                                      x[,k], z, tau, u[k], scale)}))/N
  }
  noise
}

#####
# Calculation of the latent trace in the M-step based off Listgarten et al
#     paper
#
m_step_z <- function(K,M,Q,state_prob,x,u,scale,noise,lambda,periodic){
  z <- matrix(rep(0,M),M,1)
  A <- matrix(0,M,M)
  b <- matrix(rep(0,M),M,1)
  # Adding u to penalty as suggested in Listgarten's thesis
  # to solve a degeneracy problem:
  #   Penalizing changes in z promotes small z, u can then be used
  #   to provide a shift to compensate.
  u_term <- mean(u^2)
  for (m in 1:M){
    b_term <- 0
    A_diag <- 0
    for (k in 1:K){
      for (q in 1:Q){
        state_ind <- M*(q-1)+m # Find correct row in state_prob
        not_na    <- which(!is.na(x[,k]))
        b_term    <- b_term + sum(state_prob[[k]][state_ind,not_na]*scale[q]*u[k]*x[,k][not_na]/noise)
        A_diag    <- A_diag + sum(state_prob[[k]][state_ind,not_na]*(scale[q]*u[k])^2/noise)
      }
    }
    # Penalization term half as large on ends of series due to missing neighbors
    #     Can remove this and make adjustment for periodicity when we go to
    #     periodic boundaries as we should.
    if(m==1 | m==M){
      A[m,m] <- A_diag+2*lambda*u_term
      if(periodic) A[m,m] <- A[m,m] + 2*lambda*u_term
    }else{
      A[m,m] <- A_diag+4*lambda*u_term
    }
    b[m,1] <- b_term
    if (m > 1) A[m,m-1] <- -2*lambda*u_term
    if (m < M) A[m,m+1] <- -2*lambda*u_term
  }
  if(periodic) A[1,M] <- -2*lambda*u_term
  if(periodic) A[M,1] <- -2*lambda*u_term
  
  # NOTE:
  #     When fixing periodicity, make sure to fix u-update due to u_term

  # NOTE:
  #     This solve can sometimes fail due to NA's in the state probabilities,
  #     I believe due to float underflows in the e-step. More investigation is needed.
  z <- solve(A,b)
  z
}

#####
# Optimization of scaling parameters
#
m_step_phi <- function(K,M,Q,state_prob,x,u,state_scale,z,noise,lambda,periodic){
  phi <- matrix(rep(1,Q),Q,1)
  
  for (q in 2:Q){
    num   <- 0
    denom <- 0
    for (k in 1:K){
      for (m in 1:M){
        state_ind <- M*(q-1)+m # Find correct row in state_prob
        not_na    <- which(!is.na(x[,k]))
        num   <- num + sum(state_prob[[k]][state_ind,not_na]*z[m]*u[k]*x[,k][not_na]/noise)
        denom <- denom + sum(state_prob[[k]][state_ind,not_na]*(z[m]*u[k])^2/noise)
      }
    }
    phi[q] <- num/denom
  }
  
  if ( any(is.na(phi)) ) return(phi)
  if ( any(sort(phi) != phi) || min(phi) < 1) return(state_scale)
  phi
}

expected_likelihood <- function(K,M,Q,state_prob,x,z,u,scale,noise,lambda,periodic){
  l_like <- 0
  u_term <- mean(u^2)
  for (m in 1:M){
    for (k in 1:K){
      for (q in 1:Q){
        state_ind <- M*(q-1)+m # Find correct row in state_prob
        not_na    <- which(!is.na(x[,k]))
        l_like    <- l_like - sum(state_prob[[k]][state_ind,not_na]*(log(noise)/2+(x[,k][not_na]-z[m]*scale[q]*u[k])^2/2/noise))
      }
    }
    if (m!=M) l_like <- l_like - lambda*u_term*(z[m+1]-z[m])^2
    if(periodic && m==M) l_like <- l_like - lambda*u_term*(z[1]-z[m])^2
  }
  l_like
}

#####
# Update of the uniform scaling parameter u_k using equation on Page 4 of
#     2004 paper.
#
m_step_u <- function(K,M,Q,state_prob,x,z,tau,scale,n_k,lambda,noise,periodic,training = TRUE){
  u <- rep(0,K)
  for (k in 1:K){
    num    <- 0
    denom  <- 0
    not_na <- which(!is.na(x[,k]))
    for (m in 1:M){
      for (q in 1:Q){
        state_ind <- M*(q-1)+m # Find correct row in state_prob
        num   <- num + z[m]*scale[[q]]*sum(state_prob[[k]][state_ind,not_na]*x[,k][not_na])/noise
        denom <- denom + (z[m]*scale[[q]])^2*(sum(state_prob[[k]][state_ind,not_na]))/noise
      }
      # Added term to deal with the u-dependence of the penalty
      if(m != M && training) denom <- denom + 2*lambda*(z[m+1] - z[m])^2/K
      if(periodic && m==M && training) denom <- denom + 2*lambda*(z[1] - z[M])^2/K
    }
    u[k] <- num/denom
  }
  u
}


#####
# Core function for running the alignment algorithm.
#
align_series_EM <- function(t,x,
                            J=3,
                            upsample_factor=2,
                            buffer=NULL, 
                            tol = 1.e-5,
                            transition_tau = NULL,
                            transition_scale = c(0.5,0.25),
                            noise=NULL, 
                            init_t_range = NULL,
                            scales = 1, 
                            iters = 20,
                            lambda=0,
                            init_z = NULL,
                            n_times = NULL,
                            parallel_cores=1,
                            periodic=FALSE,
                            get_viterbi=TRUE){
  # Set parameter defaults
  if(is.null(noise)) noise <- (0.15*diff(range(unlist(x))))^2
  if(is.null(buffer)) buffer <-floor(0.05*max(lengths(x)))
  if(is.null(transition_tau)) transition_tau <- rep(1/J,J)

  # Error check the number of cores
  if( parallel_cores != as.integer(parallel_cores) || parallel_cores < 1) stop('Number of cores must be a positive integer')
  if(parallel_cores > 1) {
    library(parallel)
    num_cores <- detectCores()
    if (is.na(num_cores)) {
      print("Number of detected cores == NA. Running in serial.")
      parallel_cores <- 1
    } else {
      if( parallel_cores > num_cores) parallel_cores <- num_cores
      if( parallel_cores == num_cores) print("Using the maximum number of cores")
    }
  }
  
  if (is.null(n_times)){
    tmp <- unique(unlist(t))
    n_times <- length(tmp[!is.na(tmp)])  # Length of time series assuming full set of observations
  }
  if(is.null(init_t_range)) init_t_range <- 1:max(floor(0.1*n_times),J)
  K <- length(x)               # Number of replicate time series
  M <- floor(upsample_factor*n_times)+2*buffer # Number of latent times
  N_k <- lengths(x)            # lengths of each replicate series
  Q <- length(scales)          # Number of latent scales
  ind <- which.max(lengths(x)) # Find longest series for initializing latent profile
  states <- expand.grid(t=1:M,q=1:Q) # time, scale state pairs
  
  state_tau <- 1:M             # Time indices for latent profile
  state_scale <- scales        # Possible states for scale parameter
  log_like <- c()              # Store the log-likelihood
  
  latent_tau <- seq(min(unlist(t)),max(unlist(t)),length = M) # Latent times
  if (is.null(init_z)) init_z <- x[[1]]
  latent_z <- up_sample(init_z,M)  # Latent profile vector (z)  

  uniform_scale <- rep(1,K)    # Global scaling for each replicate series (u)
  for ( ii in 1:K ) uniform_scale[ii] <- mean(x[[ii]])/mean(latent_z)

  # Create matrix for even-spacing time series with NA values for unobserved times 
  #     x_list and t_list contain raw times and observations for uneven sampling
  x_list <- x
  x <- matrix(NA,nrow=n_times,ncol=K)
  if(!is.list(t)){
    if(!all(lengths(x_list)==length(t))) stop('t is a vector but not all series in x have same length as t')
    for (i in 1:K) x[round(t),i] <- x_list[[i]]
  } else{
    for (i in 1:K) x[round(t[[i]]),i] <- x_list[[i]]
  }
  t_list <- t
  t <- 1:n_times
  # Normalize transition probabilities for tau and scale states
  transition_tau <- transition_tau/sum(transition_tau)
  transition_scale <- transition_scale/sum(transition_scale*c(1,2))

  # Build the list of state_probability matrices
  state_prob <- list(matrix(1,M*Q,n_times)/(M*Q))
  for (k in 2:K)   state_prob[[k]] <- matrix(1,M*Q,n_times)/(M*Q)  

  # Perform the E-M algorithm to fit the model parameters and obtain state probabilities
  cat("Model Fitting Progress\n Iteration:")
  old_log_like <- -Inf
  noise_check_passed <- FALSE
  for (iter in 1:iters){
    if(iter%%5 == 0) cat(" ",iter)    
    old_noise <- noise
    old_latent_z <- latent_z
    old_uniform_scale <- uniform_scale
    state_prob <- e_step(K,M,Q,x,latent_z,latent_tau,state_scale,state_prob,
                               uniform_scale,transition_tau,transition_scale,
                               noise, init_t_range,parallel_cores,periodic)
    if(!noise_check_passed){ 
      if( any(unlist(lapply(state_prob,function(x) any(apply(x,2,function(y) all(y==0))))))){
        noise <- noise*10
        next
      } else{
        noise_check_passed <- TRUE
      }
    }
    noise <- m_step_noise(K,M,Q,state_prob,x,latent_z,state_tau,
                          uniform_scale,state_scale,N_k)
    latent_z <- m_step_z(K,M,Q,state_prob,x,uniform_scale,
                         state_scale,noise,lambda,periodic)
    if (Q > 1) state_scale <- m_step_phi(K,M,Q,state_prob,x,uniform_scale,
                         state_scale,latent_z,noise,lambda,periodic)
    uniform_scale <- m_step_u(K,M,Q,state_prob,x,latent_z,state_tau,
                              state_scale,N_k,lambda,noise,periodic,
                              training = TRUE)
    change_z <- abs(latent_z-old_latent_z)
    change_noise <- abs(noise-old_noise)
    change_u <- abs(uniform_scale-old_uniform_scale)
    log_like <- c(log_like,expected_likelihood(K,M,Q,state_prob,x,latent_z,uniform_scale,state_scale,noise,lambda,periodic))
    if(any(is.nan(change_u))) stop('NaN detected after M-step. Potential float underflow. Attempt different initialization; if problem persists check model parameters.')
    if((log_like[iter] - old_log_like)/abs(log_like[iter]) < tol) break
    old_log_like <- log_like[iter]
    # if(max(c(change_noise,change_z,change_u))<tol) break    
  }
  if(iters>1) cat('\n')
  if(iter==iters) print("Failed to Converge!") # Warn user that E-M did not converge

  max_prob <- sapply(state_prob,function(x){apply(x,2,max)})
  max_prob_loc <- sapply(state_prob,function(x){apply(x,2,which.max)})

  # Get maximum probability path using the Viterbi algorithm
  if(get_viterbi){
    viterbi_path <- get_viterbi_path(K,M,Q,x,latent_z,latent_tau,state_scale,
                                     uniform_scale,transition_tau,transition_scale,
                                     noise,init_t_range,parallel_cores,periodic)
    viterbi_z <- apply(viterbi_path,2,function(x){latent_z[states$t[x]]})
    viterbi_tau <- apply(viterbi_path,2,function(x){latent_tau[states$t[x]]})
  }else{
    viterbi_path <- NULL
    viterbi_z <- NULL
    viterbi_tau <- NULL
  }
  # Return list of values
  # NOTE: May want to make this into a custom R object. May not have time
  #           before I go though. Not sure how important it is.
  list(z=latent_z,tau=latent_tau,u=uniform_scale,noise=noise, scales=state_scale, states=states, max_prob=max_prob, max_prob_loc = max_prob_loc,
       viterbi_path=viterbi_path,viterbi_z=viterbi_z,viterbi_tau=viterbi_tau,lambda=lambda,log_like=log_like)
}

#####
# Run a restricted EM to obtain the predicted alignment
#   with the latent profile and noise parameter fixed
#
# NOTE: J should not need to be settable here. Should probably
#           include it from the reference alignment model
predict_series_EM <- function(t,x,aligned_obj,
                              tol = 1.e-5,
                              J=3,
                              init_t_range=NULL,
                              transition_tau = NULL,
                              transition_scale = c(0.5,0.25),
                              iters = 20, 
                              parallel_cores=1, 
                              periodic=FALSE,
                              get_viterbi=TRUE){
  # Get variables from previous alignment
  z <- aligned_obj$z
  noise <- aligned_obj$noise
  scales <- aligned_obj$scales
  n_times <- nrow(aligned_obj$viterbi_z)
  lambda <- aligned_obj$lambda
  
  # Error check the number of cores
  if(parallel_cores != as.integer(parallel_cores) || parallel_cores < 1) stop('Number of cores must be a positive integer')
  if(parallel_cores > 1) {
    library(parallel)
    num_cores <- detectCores()
    if (is.na(num_cores)) {
      print("Number of detected cores == NA. Running in serial.")
      parallel_cores <- 1
    } else {
      if( parallel_cores > num_cores) parallel_cores <- num_cores
      if( parallel_cores == num_cores) print("Using the maximum number of cores")
    }
  }

  if(is.null(init_t_range)) init_t_range <- 1:max(floor(0.1*n_times),J)
  K <- length(x)               # Number of replicate time series
  M <- length(z)               # Number of latent times
  N_k <- lengths(x)            # lengths of each replicate series
  Q <- length(scales)          # Number of latent scales
  
  if(!is.list(x)) stop('Function requires at least two replicate series for prediction')
  
  state_tau <- 1:M             # Time indices for latent profile
  state_scale <- scales        # Possible states for scale parameter
  if(is.null(transition_tau)) transition_tau <- rep(1/J,J)
  log_like <- c()
  
  latent_z <- z                # Latent profile vector (z)
  latent_tau <- seq(min(unlist(t)),max(unlist(t)),length = M) # Latent times

  uniform_scale <- rep(1,K)    # Global scaling for each replicate series (u)
  for ( ii in 1:K ) uniform_scale[ii] <- mean(x[[ii]])/mean(latent_z)

  # Create matrix for even-spacing time series with NA values for unobserved times 
  #     x_list and t_list contain raw times and observations for uneven sampling
  x_list <- x
  x <- matrix(NA,nrow=n_times,ncol=K)
  if(!is.list(t)){
    if(!all(lengths(x_list)==length(t))) stop('t is a vector but not all series in x have same length as t')
    for (i in 1:K){x[t,i] <- x_list[[i]]}
  } else{
    for (i in 1:K){x[t[[i]],i] <- x_list[[i]]}    
  }
  t <- 1:n_times
  
  # Normalize transition probabilities for tau and scale states
  transition_tau <- transition_tau/sum(transition_tau)
  transition_scale <- transition_scale/sum(transition_scale*c(1,2))
  
  # Build the list of state_probability matrices
  state_prob <- list(matrix(1,M*Q,n_times)/(M*Q)) 
  for (k in 2:K)   state_prob[[k]] <- matrix(1,M*Q,n_times)/(M*Q)  
  
  # Perform the restricted E-M algorithm to fit the series specific model parameters 
  #     and obtain state probabilities
  cat("Model Fitting Progress\n Iteration:")
  old_log_like <- -Inf
  for (iter in 1:iters){
    if(iter%%5 == 0) cat(" ",iter)
    old_uniform_scale <- uniform_scale
    mult <- ifelse(iter==1,10,1)

    state_prob <- e_step(K,M,Q,x,latent_z,latent_tau,state_scale,state_prob,
                         uniform_scale,transition_tau,transition_scale,
                         mult*noise,init_t_range,parallel_cores,periodic)
    uniform_scale <- m_step_u(K,M,Q,state_prob,x,latent_z,state_tau,
                              state_scale,N_k,lambda,mult*noise,periodic,
                              training = FALSE)
    change_u <- abs(uniform_scale-old_uniform_scale)
    log_like <- c(log_like,expected_likelihood(K,M,Q,state_prob,x,latent_z,uniform_scale,state_scale,mult*noise,lambda,periodic))
    if(any(is.nan(change_u))) stop('NaN detected after M-step')
    if((log_like[iter] - old_log_like)/abs(log_like[iter]) < tol) break
    old_log_like <- log_like[iter]
    # if(max(change_u)<1.e-5) break
  }
  if(iters>1) cat('\n')
  if(iter==iters) print("Failed to Converge!") # Warn user that E-M did not converge
  
  # Get maximum probability path using the Viterbi algorithm
  if(get_viterbi){
    states <- expand.grid(t=1:M,q=1:Q)
    viterbi_path <- get_viterbi_path(K,M,Q,x,latent_z,latent_tau,state_scale,
                                     uniform_scale,transition_tau,transition_scale,
                                     noise,init_t_range,parallel_cores,periodic)
    viterbi_z <- apply(viterbi_path,2,function(x){latent_z[states$t[x]]})
    viterbi_tau <- apply(viterbi_path,2,function(x){latent_tau[states$t[x]]})
  }else{
    viterbi_path <- NULL
    viterbi_z <- NULL
    viterbi_tau <- NULL
  }
  # Return list of values
  # NOTE: May want to make this into a custom R object. May not have time
  #           before I go though. Not sure how important it is.
  list(z=latent_z,tau=latent_tau,u=uniform_scale,noise=noise, scales=scales, states=states,
       viterbi_path=viterbi_path,viterbi_z=viterbi_z,viterbi_tau=viterbi_tau,log_like=log_like)
}

get_residuals <- function(t,x,aligned_obj, 
                          standardize=FALSE){
  n_series <- length(aligned_obj$u)
  fit_z <- aligned_obj$viterbi_z
  fit_t <- aligned_obj$viterbi_tau
  if(is.list(t)){
    model_fit <- lapply(1:n_series,function(i){fit_z[t[[i]],i]*aligned_obj$u[i]*
                                                 aligned_obj$scales[aligned_obj$states$q[aligned_obj$viterbi_path[t[[i]],i]]]})
  }else {
    model_fit <- lapply(1:n_series,function(i){fit_z[t,i]*aligned_obj$u[i]*
                                                 aligned_obj$scales[aligned_obj$states$q[aligned_obj$viterbi_path[t,i]]]})
  }
  resid <- mapply(function(x,y){x-y}, x,model_fit,SIMPLIFY = F)
  if(standardize) resid <- lapply(resid,function(x){(x-mean(unlist(resid)))/sqrt(aligned_obj$noise)})
  resid
}

#####
# Run a restricted EM to obtain the predicted alignment
#   with the latent profile and noise parameter fixed
#
# NOTE: J should not need to be settable here. Should probably
#           include it from the reference alignment model
online_align <- function(t,x,aligned_obj,
                         tol = 1.e-5,
                         J   = 3,
                         init_t_range   = NULL,
                         transition_tau = NULL,
                         transition_scale = c(0.5,0.25),
                         iters = 20, 
                         parallel_cores = 1, 
                         periodic = FALSE,
                         get_viterbi = TRUE){
  # Get variables from previous alignment
  z       <- aligned_obj$z
  noise   <- aligned_obj$noise
  scales  <- aligned_obj$scales
  n_times <- nrow(aligned_obj$viterbi_z)
  lambda  <- aligned_obj$lambda
  
  # Error check the number of cores
  if( parallel_cores != as.integer(parallel_cores) || parallel_cores < 1) stop('Number of cores must be a positive integer')
  if(parallel_cores > 1) {
    library(parallel)
    num_cores <- detectCores()
    if (is.na(num_cores)) {
      print("Number of detected cores == NA. Running in serial.")
      parallel_cores <- 1
    } else {
      if( parallel_cores > num_cores) parallel_cores <- num_cores
      if( parallel_cores == num_cores) print("Using the maximum number of cores")
    }
  }

  if(is.null(init_t_range)) init_t_range <- 1:max(floor(0.1*n_times),J)
  M   <- length(z)             # Number of latent times
  N_k <- length(x)             # lengths of each replicate series
  Q   <- length(scales)        # Number of latent scales
  
  state_tau   <- 1:M           # Time indices for latent profile
  state_scale <- scales        # Possible states for scale parameter
  if(is.null(transition_tau)) transition_tau <- rep(1/J,J)
  log_like <- c()
  
  uniform_scale <- 1           # Global scaling for each replicate series (u)
  latent_z      <- z           # Latent profile vector (z)
  latent_tau    <- seq(min(unlist(t)),max(unlist(t)),length = M) # Latent times

  # Create matrix for even-spacing time series with NA values for unobserved times 
  #     x_list and t_list contain raw times and observations for uneven sampling
  x_list <- x
  x      <- matrix(NA,nrow=n_times,ncol=1)
  x[t,1] <- x_list
  t      <- 1:n_times
  
  # Normalize transition probabilities for tau and scale states
  transition_tau   <- transition_tau/sum(transition_tau)
  transition_scale <- transition_scale/sum(transition_scale*c(1,2))
  
  # Build the list of state_probability matrices
  state_prob <- list(matrix(1,M*Q,n_times)/(M*Q)) 
  
  # Perform the restricted E-M algorithm to fit the series specific model parameters 
  #     and obtain state probabilities
  cat("Model Fitting Progress\n Iteration:")
  old_log_like <- -Inf
  for (iter in 1:iters){
    if(iter%%5 == 0) cat(" ",iter)
    old_uniform_scale <- uniform_scale

    state_prob <- e_step(K=1,M,Q,x,latent_z,latent_tau,state_scale,state_prob,
                         uniform_scale,transition_tau,transition_scale,
                         noise,init_t_range,parallel_cores,periodic)
    uniform_scale <- m_step_u(K=1,M,Q,state_prob,x,latent_z,state_tau,
                              state_scale,N_k,lambda,noise,periodic,
                              training=FALSE)
    change_u <- abs(uniform_scale-old_uniform_scale)
    log_like <- c(log_like,expected_likelihood(K=1,M,Q,state_prob,x,latent_z,uniform_scale,state_scale,noise,lambda,periodic))
    if(any(is.nan(change_u))) stop('NaN detected after M-step')
    if((log_like[iter] - old_log_like)/abs(log_like[iter]) < tol) break
    old_log_like <- log_like[iter]
  }
  if(iters>1) cat('\n')
  if(iter==iters) print("Failed to Converge!") # Warn user that E-M did not converge
  
  # Get maximum probability path using the Viterbi algorithm
  if(get_viterbi){
    states       <- expand.grid(t=1:M,q=1:Q)
    viterbi_path <- get_viterbi_path(K=1,M,Q,x,latent_z,latent_tau,state_scale,
                                     uniform_scale,transition_tau,transition_scale,
                                     noise,init_t_range,parallel_cores,periodic)
    viterbi_z   <- apply(viterbi_path,2,function(x){latent_z[states$t[x]]})
    viterbi_tau <- apply(viterbi_path,2,function(x){latent_tau[states$t[x]]})
  }else{
    viterbi_path <- NULL
    viterbi_z    <- NULL
    viterbi_tau  <- NULL
  }
  # Return list of values
  # NOTE: May want to make this into a custom R object. May not have time
  #           before I go though. Not sure how important it is.
  list(z=latent_z,tau=latent_tau,u=uniform_scale,noise=noise, scales=scales, states=states,
       viterbi_path=viterbi_path,viterbi_z=viterbi_z,viterbi_tau=viterbi_tau,log_like=log_like)
}

online_residuals <- function(t,x,aligned_obj, 
                             standardize=FALSE){
  fit_z     <- aligned_obj$viterbi_z
  model_fit <- fit_z[t,1]*aligned_obj$u[1]*
                          aligned_obj$scales[aligned_obj$states$q[aligned_obj$viterbi_path[t,1]]]
  resid <- x-model_fit
  if(standardize) resid <- (resid-mean(resid))/sqrt(aligned_obj$noise)
  resid
}

init_filter_state <- function(t,x,aligned_obj){
  tau   <- aligned_obj$states[aligned_obj$viterbi_path[t,],1]
  scale <- aligned_obj$states[aligned_obj$viterbi_path[t,],2]
  resid <- online_residuals(t,x,aligned_obj)
  u     <- aligned_obj$u

  data.frame(tau=tau,scale=scale,residual=resid,u=u,t=t,x=x)
}

filter_align <- function(t,x,aligned_obj, filter_state,
                            J=3,
                            transition_tau   = NULL,
                            transition_scale = c(0.5,0.25), 
                            periodic=FALSE){
  z       <- aligned_obj$z
  noise   <- aligned_obj$noise
  scales  <- aligned_obj$scales
  states  <- aligned_obj$states
  n_obs   <- length(filter_state$tau)
  u       <- filter_state$u[1]
  M       <- length(z)
  Q       <- length(scales)

  if(is.null(transition_tau)) transition_tau <- rep(1/J,J)
  if(Q>1){
    transition_scale <- transition_scale/sum(transition_scale*c(1,2))
    transition_scale <- c(transition_scale[2],transition_scale)
  }

  last_tau   <- filter_state$tau[n_obs] 
  last_scale <- filter_state$scale[n_obs]
  
  p_tau   <- 1.
  p_scale <- rep(0,Q)

  p_scale[last_scale] <- 1.
  if(Q>1){
    t_scale_matrix      <- diag(Q)*transition_scale[2]

    t_scale_matrix[matrix(c(2:Q,1:(Q-1)),Q-1,2)] <- transition_scale[1]
    t_scale_matrix[matrix(c(1:(Q-1),2:Q),Q-1,2)] <- transition_scale[1]
    t_scale_matrix <- t(apply(t_scale_matrix,1,function(x) x/sum(x)))
  } else{
    t_scale_matrix <- 1
  }

  possible_tau   <- last_tau
  possible_scale <- 1:Q
  
  for (steps in 1:(t-filter_state$t[n_obs])){

    possible_tau   <- (possible_tau[1]+1):(possible_tau[length(possible_tau)]+J)
    p_tmp <- p_tau
    p_tau <- rep(0,length(possible_tau))
    for (val in 1:length(p_tmp)) {
      p_tau[val+0:(J-1)] <- p_tau[val+0:(J-1)] + p_tmp[val]*transition_tau
    }

    p_scale <- t(t_scale_matrix) %*% p_scale

  }

  if(periodic){
    possible_tau <- (possible_tau-1) %% M + 1
  }else{
    if(possible_tau[1] > M) stop('No states to transition into')
    tmp <- which(possible_tau > M)
    possible_tau   <- possible_tau[-tmp]
    transition_tau <- transition_tau[-tmp]
    transition_tau <- transition_tau/sum(transition_tau)
  }

  possible_states  <- expand.grid(possible_tau,possible_scale)
  transition_prob  <- expand.grid(p_tau,p_scale)
  transition_prob  <- apply(transition_prob,1,prod)
  state_mean       <- z[possible_states[,1]]*u*scales[possible_states[,2]]
  state_likelihood <- dnorm(x,state_mean,sqrt(noise))
  state_prob       <- state_likelihood*transition_prob

  best_ind   <- which.max(state_prob)
  best_state <- possible_states[best_ind,]
  best_scale <- possible_states[best_ind,2]
  best_tau   <- possible_states[best_ind,1]
  best_mean  <- state_mean[best_ind]

  # new_state <- which(apply(states,1,check_states,best_state))
  residual  <- x - best_mean
  
#if(abs(residual) > 1.0){
#  print(state_prob)
#  print(possible_states)
#  print(state_mean)
#  print(t)
#  print(x)
#  print(c(z))
#}
  filter_state[n_obs+1,] <- c(best_tau,best_scale,residual,u,t,x)

  filter_state
}
 
init_particle_filter <- function(n_particles,n_obs,init_tau,init_scale,
                                 init_p_tau   = NULL,
                                 init_p_scale = NULL,
                                 u_range = c(0.97,1.03)){
  n_tau   <- length(init_tau)
  n_scale <- length(init_scale)

  if (is.null(init_p_tau))   init_p_tau   <- rep(1/n_tau, n_tau)
  if (is.null(init_p_scale)) init_p_scale <- rep(1/n_scale,n_scale)

  taus   <- sample(init_tau,   n_particles, replace=TRUE, prob=init_p_tau)
  scales <- sample(init_scale, n_particles, replace=TRUE, prob=init_p_scale)
  us     <- runif(n_particles, u_range[1],  u_range[2])

  particles <- list() 
  for ( ii in 1:n_particles){
    particles[[ii]] <- data.frame(tau      = rep(taus[ii],   n_obs+1),
                                  scale    = rep(scales[ii], n_obs+1),
                                  u        = rep(us[ii],     n_obs+1),
                                  t        = rep(0,          n_obs+1),
                                  x        = rep(0,          n_obs+1),
                                  residual = rep(0.,         n_obs+1))
  }

  particles
}

particle_align <- function(t,x,aligned_obj, particles,ind,
                            J=3,
                            transition_tau   = NULL,
                            transition_scale = c(0.5,0.25), 
                            periodic=FALSE){
  z           <- aligned_obj$z
  noise       <- aligned_obj$noise
  scales      <- aligned_obj$scales
  states      <- aligned_obj$states
  M           <- length(z)
  Q           <- length(scales)
  n_particles <- length(particles)
  weights     <- rep(1/n_particles,n_particles)

  if(is.null(transition_tau)) transition_tau <- rep(1/J,J)
  if(Q>1){
    transition_scale <- transition_scale/sum(transition_scale*c(1,2))
    transition_scale <- c(transition_scale[2],transition_scale)
  }

 

  if(Q>1){
    t_scale_matrix <- diag(Q)*transition_scale[2]

    t_scale_matrix[matrix(c(2:Q,1:(Q-1)),Q-1,2)] <- transition_scale[1]
    t_scale_matrix[matrix(c(1:(Q-1),2:Q),Q-1,2)] <- transition_scale[1]

    t_scale_matrix <- t(apply(t_scale_matrix,1,function(x) x/sum(x)))
  } else{
    t_scale_matrix <- 1
  }

  possible_scale <- 1:Q
  
  for (ii in 1:n_particles){
    p_tau        <- 1.
    last_tau     <- particles[[ii]]$tau[ind-1] 
    last_scale   <- particles[[ii]]$scale[ind-1]
    possible_tau <- last_tau
    p_scale      <- rep(0,Q)
    p_scale[last_scale] <- 1.

    for (steps in 1:(t-particles[[ii]]$t[ind-1])){
     
      u <- particles[[ii]]$u[1]

      possible_tau   <- (possible_tau[1]+1):(possible_tau[length(possible_tau)]+J)

      p_tmp <- p_tau
      p_tau <- rep(0,length(possible_tau))

      for (val in 1:length(p_tmp)) {
        p_tau[val+0:(J-1)] <- p_tau[val+0:(J-1)] + p_tmp[val]*transition_tau
      }

      p_scale <- t(t_scale_matrix) %*% p_scale
    }

    if(periodic){
      possible_tau <- (possible_tau-1) %% M + 1
    }else{
      if(possible_tau[1] > M) stop('No states to transition into')
      tmp <- which(possible_tau > M)
      possible_tau   <- possible_tau[-tmp]
      transition_tau <- transition_tau[-tmp]
      transition_tau <- transition_tau/sum(transition_tau)
    }

    possible_states  <- expand.grid(possible_tau,possible_scale)
    transition_prob  <- expand.grid(p_tau,p_scale)
    transition_prob  <- apply(transition_prob,1,prod)
    state_mean       <- z[possible_states[,1]]*u*scales[possible_states[,2]]
    log_likelihood   <- dnorm(x,state_mean,sqrt(noise),log=TRUE)
    
    new_ind     <- sample(1:nrow(possible_states),1,prob=transition_prob)
    new_state   <- possible_states[new_ind,]
    new_scale   <- possible_states[new_ind,2]
    new_tau     <- possible_states[new_ind,1]
    new_mean    <- state_mean[new_ind]
    residual    <- x - new_mean
    weights[ii] <- log_likelihood[new_ind]

    particles[[ii]][ind,] <- c(new_tau, new_scale,u,t,x,residual)
  } 

  weights <- exp(weights - max(weights))
  if (all(weights==0)) print("Warning!! particle weights all zero")
  resample_inds <- sample(1:n_particles, n_particles,replace=TRUE, prob=weights+.Machine$double.xmin)
  particles     <- particles[resample_inds]

  particles
}
 
