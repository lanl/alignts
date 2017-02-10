######
# Wrappers for FORTRAN functions to speed up
#   the code for aligning replicate time series
#   and inferring the latent series profile 
#   based on Listgarten, Neal, Roweis, Emili
#
# Required libraries:
#   NONE
#
######

#####
# Wrapper for FORTRAN function for calculating E-step. Contributes huge speed-up
#     necessary to use algorithm on full set of data. Function is also 
#     parallelized using OpenMP to contribute further performance boosts.
#
e_step <- function(K,M,Q,x,z,tau,scale,state_prob,u,transition_tau,
                   transition_scale,noise,init_t_range,parallel_cores,periodic){
  if (parallel_cores > 1){
    if(periodic){
      e_step_call <- e_step_periodic
    } else{
      e_step_call <- e_step_parallel
    }
  }else{
    if(periodic){ 
      e_step_call <- e_step_serial_periodic
    } else{
      e_step_call <- e_step_serial
    }
  }
  states <- expand.grid(t=1:M,q=1:Q)
  s <- nrow(states)
  N <- nrow(x)
  init_n <- length(init_t_range)
  init_t_range <- ((init_t_range - 1) %% M) + 1
  # R-FORTRAN interface does not like passing NA, so passing bad value in place
  #     unlikely to ever come across legit data that small.
  x[which(is.na(x))] <- -1.e25
  n_tau <- length(transition_tau)
  n_scale <- length(transition_scale)
  state_prob <- rep(0,s*N*K)
  retdata <- .Fortran(e_step_call,
                      ncores = as.integer(parallel_cores),
                      K=as.integer(K),
                      M=as.integer(M),
                      Q=as.integer(Q),
                      N=as.integer(N),
                      n_tau=as.integer(n_tau),
                      n_scale=as.integer(n_scale),
                      x = as.matrix(x),
                      z = as.double(z),
                      states=as.matrix(states),
                      noise = as.double(noise),
                      init_n = as.integer(init_n),
                      init_t_range = as.integer(init_t_range),
                      u=as.double(u),
                      phi=as.double(scale),
                      t_tau=as.double(transition_tau),
                      t_scale=as.double(transition_scale),
                      forward_mat=as.double(rep(0,(M*Q+1)*(N+1)*K)),
                      backward_mat=as.double(rep(0,(M*Q+1)*(N+1)*K)),
                      transition_prob=as.double(rep(0,M*Q*n_tau*3)),
                      transition_inds=as.integer(rep(0,M*Q*n_tau*3*2)),
                      state_prob=as.double(state_prob))$state_prob
  # Reshape output from FORTRAN routine into the necessary list of matrices
  state_prob <- list()
  for (i in 1:K){state_prob[[i]] <- matrix(retdata[1:(s*N)+(i-1)*s*N],s,N)}
  state_prob
}

#####
# Wrapper for FORTRAN function for performing Viterbi algorithm. Contributes 
#     huge speed-up necessary to use algorithm on full set of data. Function is
#     also parallelized using OpenMP to contribute further performance boosts.
#
get_viterbi_path <- function(K,M,Q,x,z,tau,scale,u,transition_tau,
                   transition_scale,noise,init_t_range,parallel_cores,periodic){
  if(parallel_cores){
    if(periodic){
      get_viterbi_path_call <- get_viterbi_path_periodic
    } else{
      get_viterbi_path_call <- get_viterbi_path_parallel
    }
  } else {
    if(periodic){
      get_viterbi_path_call <- get_viterbi_path_serial_periodic
    } else{
      get_viterbi_path_call <- get_viterbi_path_serial
    }
  }
  states <- expand.grid(t=1:M,q=1:Q)
  s <- nrow(states)
  N <- nrow(x)
  init_n <- length(init_t_range)
  init_t_range <- ((init_t_range - 1) %% M) + 1
  # R-FORTRAN interface does not like passing NA, so passing bad value in place
  #     unlikely to ever come across legit data that small.
  x[which(is.na(x))] <- -1.e25
  n_tau <- length(transition_tau)
  n_scale <- length(transition_scale)
  best_path <- rep(0,N*K)
  retdata <- .Fortran(get_viterbi_path_call,
                      ncores = as.integer(parallel_cores),
                      K=as.integer(K),
                      M=as.integer(M),
                      Q=as.integer(Q),
                      N=as.integer(N),
                      n_tau=as.integer(n_tau),
                      n_scale=as.integer(n_scale),
                      x = as.matrix(x),
                      z = as.double(z),
                      states = as.matrix(states),
                      noise  = as.double(noise),
                      init_n = as.integer(init_n),
                      init_t_range = as.integer(init_t_range),
                      u=as.double(u),
                      phi=as.double(scale),
                      t_tau=as.double(transition_tau),
                      t_scale=as.double(transition_scale),
                      V_mat=as.double(rep(0,M*Q*N*K)),
                      possible_paths=as.integer(rep(0,M*Q*N*K)),
                      best_path=as.integer(best_path))
  # Reshape output from FORTRAN routine into the necessary matrix
  best_path <- matrix(retdata$best_path,N,K)
  best_path
}
