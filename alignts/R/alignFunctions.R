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
    tmp <- x[sort(inds)]
  }
  tmp
}

#####
# Calculated the probability-weighted SSE given a state and set
#     of state probabilities
#
state_sse <- function(q, t, state_prob, x, z, tau, u, scale){
  not_na <- which(!is.na(x))
  sum(state_prob[not_na]*(x[not_na]-z[tau[t]]*u*scale[q])^2)
}

#####
# Update sigma^2 using equation on Page 4 of paper
#
m_step_noise <- function(K,M,Q,state_prob,x,z,tau,u,scale,n_k){
  states <- expand.grid(t=1:M,q=1:Q)
  noise <- 0
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
        not_na <- which(!is.na(x[,k]))
        b_term <- b_term + sum(state_prob[[k]][state_ind,not_na]*scale[q]*u[k]*x[,k][not_na]/noise)
        A_diag <- A_diag + sum(state_prob[[k]][state_ind,not_na]*(scale[q]*u[k])^2/noise)
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
m_step_phi <- function(K,M,Q,state_prob,x,u,z,noise,lambda,periodic){
  phi <- matrix(rep(1,Q),Q,1)
  
  for (q in 2:Q){
    num <- 0
    denom <- 0
    for (k in 1:K){
      for (m in 1:M){
        state_ind <- M*(q-1)+m # Find correct row in state_prob
        not_na <- which(!is.na(x[,k]))
        num <- num + sum(state_prob[[k]][state_ind,not_na]*z[m]*u[k]*x[,k][not_na]/noise)
        denom <- denom + sum(state_prob[[k]][state_ind,not_na]*(z[m]*u[k])^2/noise)
      }
    }
    phi[q] <- num/denom
  }
  
  phi
}

expected_likelihood <- function(K,M,Q,state_prob,x,z,u,scale,noise,lambda,periodic){
  l_like <- 0
  u_term <- mean(u^2)
  for (m in 1:M){
    for (k in 1:K){
      for (q in 1:Q){
        state_ind <- M*(q-1)+m # Find correct row in state_prob
        not_na <- which(!is.na(x[,k]))
        l_like <- l_like - sum(state_prob[[k]][state_ind,not_na]*(log(noise)/2+(x[,k][not_na]-z[m]*scale[q]*u[k])^2/2/noise))
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
m_step_u <- function(K,M,Q,state_prob,x,z,tau,scale,n_k,lambda,noise,periodic){
  u <- rep(0,K)
  for (k in 1:K){
    num <- 0
    denom <- 0
    not_na <- which(!is.na(x[,k]))
    for (m in 1:M){
      for (q in 1:Q){
        state_ind <- M*(q-1)+m # Find correct row in state_prob
        num <- num + z[m]*scale[[q]]*sum(state_prob[[k]][state_ind,not_na]*x[,k][not_na])/noise
        denom <- denom + (z[m]*scale[[q]])^2*(sum(state_prob[[k]][state_ind,not_na]))/noise
      }
      # Added term to deal with the u-dependence of the penalty
      if(m != M) denom <- denom + 2*lambda*(z[m+1] - z[m])^2/K
      if(periodic && m==M) denom <- denom + 2*lambda*(z[1] - z[M])^2/K
    }
    u[k] <- num/denom
  }
  u
}

#####
# Core function for running the alignment algorithm.
#
align_series_EM <- function(t,x,J=5,upsample_factor=2,
                            buffer=NULL, tol = 1.e-5,
                            transition_tau = NULL,
                            transition_scale = c(0.5,0.25),
                            noise=NULL, init_t_range = NULL,
                            scales = 1, iters = 20,lambda=0,
                            init_z = NULL,
                            n_times = NULL,
                            parallel=TRUE,
                            periodic=FALSE,
                            get_viterbi=TRUE){
  # Set parameter defaults
  if(is.null(noise)) noise <- (0.15*diff(range(unlist(x))))^2
  if(is.null(buffer)) buffer <-floor(0.05*max(lengths(x)))
  if(is.null(transition_tau)) transition_tau <- rep(1/J,J)
  
  if (is.null(n_times)){
    tmp <- unique(unlist(t))
    n_times <- length(tmp[!is.na(tmp)])  # Length of time series assuming full set of observations
  }
  if(is.null(init_t_range)) init_t_range <- max(floor(0.1*n_times),J)
  K <- length(x)               # Number of replicate time series
  M <- floor(upsample_factor*n_times)+2*buffer # Number of latent times
  N_k <- lengths(x)            # lengths of each replicate series
  Q <- length(scales)          # Number of latent scales
  ind <- which.max(lengths(x)) # Find longest series for initializing latent profile
  states <- expand.grid(t=1:M,q=1:Q) # time, scale state pairs
  
  state_tau <- 1:M             # Time indices for latent profile
  state_scale <- scales        # Possible states for scale parameter
  log_like <- c()              # Store the log-likelihood
  
  uniform_scale <- rep(1,K)    # Global scaling for each replicate series (u)
  latent_tau <- seq(min(unlist(t)),max(unlist(t)),length = M) # Latent times
  if (is.null(init_z)) init_z <- x[[1]]
  latent_z <- up_sample(init_z,M)  # Latent profile vector (z)  
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
  for (iter in 1:iters){
    if(iter%%max(1,floor(iters/10)) == 0) cat(" ",iter)    
    old_noise <- noise
    old_latent_z <- latent_z
    old_uniform_scale <- uniform_scale
    state_prob <- e_step(K,M,Q,x,latent_z,latent_tau,state_scale,state_prob,
                               uniform_scale,transition_tau,transition_scale,
                               noise, init_t_range,parallel,periodic)
    noise <- m_step_noise(K,M,Q,state_prob,x,latent_z,state_tau,
                          uniform_scale,state_scale,N_k)
    latent_z <- m_step_z(K,M,Q,state_prob,x,uniform_scale,
                         state_scale,noise,lambda,periodic)
    if (Q > 1) state_scale <- m_step_phi(K,M,Q,state_prob,x,uniform_scale,
                         latent_z,noise,lambda,periodic)
    uniform_scale <- m_step_u(K,M,Q,state_prob,x,latent_z,state_tau,
                              state_scale,N_k,lambda,noise,periodic)
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
                                     noise,init_t_range,parallel,periodic)
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
predict_series_EM <- function(t,x,aligned_obj, tol = 1.e-5,J=5,
                            init_t_range=NULL,
                            transition_tau = NULL,
                            transition_scale = c(0.5,0.25),
                            iters = 20, parallel=TRUE, 
                            periodic=FALSE,
                            get_viterbi=TRUE){
  # Get variables from previous alignment
  z <- aligned_obj$z
  noise <- aligned_obj$noise
  scales <- aligned_obj$scales
  n_times <- nrow(aligned_obj$viterbi_z)
  lambda <- aligned_obj$lambda
  
  if(is.null(init_t_range)) init_t_range <- max(floor(0.1*n_times),J)
  K <- length(x)               # Number of replicate time series
  M <- length(z)               # Number of latent times
  N_k <- lengths(x)            # lengths of each replicate series
  Q <- length(scales)          # Number of latent scales
  
  if(!is.list(x)) stop('Function requires at least two replicate series for prediction')
  
  state_tau <- 1:M             # Time indices for latent profile
  state_scale <- scales        # Possible states for scale parameter
  if(is.null(transition_tau)) transition_tau <- rep(1/J,J)
  log_like <- c()
  
  uniform_scale <- rep(1,K)    # Global scaling for each replicate series (u)
  latent_z <- z                # Latent profile vector (z)
  latent_tau <- seq(min(unlist(t)),max(unlist(t)),length = M) # Latent times

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
  old_log_like <- -Inf
  for (iter in 1:iters){
    old_uniform_scale <- uniform_scale

    state_prob <- e_step(K,M,Q,x,latent_z,latent_tau,state_scale,state_prob,
                         uniform_scale,transition_tau,transition_scale,
                         noise,init_t_range,parallel,periodic)
    uniform_scale <- m_step_u(K,M,Q,state_prob,x,latent_z,state_tau,
                              state_scale,N_k,lambda,noise,periodic)
    
    change_u <- abs(uniform_scale-old_uniform_scale)
    log_like <- c(log_like,expected_likelihood(K,M,Q,state_prob,x,latent_z,uniform_scale,state_scale,noise,lambda,periodic))
    if(any(is.nan(change_u))) stop('NaN detected after M-step')
    if((log_like[iter] - old_log_like)/abs(log_like[iter]) < tol) break
    old_log_like <- log_like[iter]
    # if(max(change_u)<1.e-5) break
  }
  if(iter==iters) print("Failed to Converge!") # Warn user that E-M did not converge
  
  # Get maximum probability path using the Viterbi algorithm
  if(get_viterbi){
    states <- expand.grid(t=1:M,q=1:Q)
    viterbi_path <- get_viterbi_path(K,M,Q,x,latent_z,latent_tau,state_scale,
                                     uniform_scale,transition_tau,transition_scale,
                                     noise,init_t_range,parallel,periodic)
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

get_residuals <- function(t,x,aligned_obj, standardize=FALSE){
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
online_align <- function(t,x,aligned_obj, tol = 1.e-5,J=5,
                            init_t_range=NULL,
                            transition_tau = NULL,
                            transition_scale = c(0.5,0.25),
                            iters = 20, parallel=TRUE, 
                            periodic=FALSE,
                            get_viterbi=TRUE){
  # Get variables from previous alignment
  z <- aligned_obj$z
  noise <- aligned_obj$noise
  scales <- aligned_obj$scales
  n_times <- nrow(aligned_obj$viterbi_z)
  lambda <- aligned_obj$lambda
  
  if(is.null(init_t_range)) init_t_range <- max(floor(0.1*n_times),J)
  M <- length(z)               # Number of latent times
  N_k <- length(x)            # lengths of each replicate series
  Q <- length(scales)          # Number of latent scales
  
  state_tau <- 1:M             # Time indices for latent profile
  state_scale <- scales        # Possible states for scale parameter
  if(is.null(transition_tau)) transition_tau <- rep(1/J,J)
  log_like <- c()
  
  uniform_scale <- 1           # Global scaling for each replicate series (u)
  latent_z <- z                # Latent profile vector (z)
  latent_tau <- seq(min(unlist(t)),max(unlist(t)),length = M) # Latent times

  # Create matrix for even-spacing time series with NA values for unobserved times 
  #     x_list and t_list contain raw times and observations for uneven sampling
  x_list <- x
  x <- matrix(NA,nrow=n_times,ncol=1)
  x[t,1] <- x_list
  t <- 1:n_times
  
  # Normalize transition probabilities for tau and scale states
  transition_tau <- transition_tau/sum(transition_tau)
  transition_scale <- transition_scale/sum(transition_scale*c(1,2))
  
  # Build the list of state_probability matrices
  state_prob <- list(matrix(1,M*Q,n_times)/(M*Q)) 
  
  # Perform the restricted E-M algorithm to fit the series specific model parameters 
  #     and obtain state probabilities
  old_log_like <- -Inf
  for (iter in 1:iters){
    old_uniform_scale <- uniform_scale

    state_prob <- e_step(K=1,M,Q,x,latent_z,latent_tau,state_scale,state_prob,
                         uniform_scale,transition_tau,transition_scale,
                         noise,init_t_range,parallel,periodic)
    uniform_scale <- m_step_u(K=1,M,Q,state_prob,x,latent_z,state_tau,
                              state_scale,N_k,lambda,noise,periodic)
    
    change_u <- abs(uniform_scale-old_uniform_scale)
    log_like <- c(log_like,expected_likelihood(K=1,M,Q,state_prob,x,latent_z,uniform_scale,state_scale,noise,lambda,periodic))
    if(any(is.nan(change_u))) stop('NaN detected after M-step')
    if((log_like[iter] - old_log_like)/abs(log_like[iter]) < tol) break
    old_log_like <- log_like[iter]
    # if(max(change_u)<1.e-5) break
  }
  if(iter==iters) print("Failed to Converge!") # Warn user that E-M did not converge
  
  # Get maximum probability path using the Viterbi algorithm
  if(get_viterbi){
    states <- expand.grid(t=1:M,q=1:Q)
    viterbi_path <- get_viterbi_path(K=1,M,Q,x,latent_z,latent_tau,state_scale,
                                     uniform_scale,transition_tau,transition_scale,
                                     noise,init_t_range,parallel,periodic)
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

online_residuals <- function(t,x,aligned_obj, standardize=FALSE){
  fit_z <- aligned_obj$viterbi_z
  model_fit <- fit_z[t,1]*aligned_obj$u[1]*
                          aligned_obj$scales[aligned_obj$states$q[aligned_obj$viterbi_path[t,1]]]
  resid <- x-model_fit
  if(standardize) resid <- (resid-mean(resid))/sqrt(aligned_obj$noise)
  resid
}


