!------------------------------------------------
! Function for calculating the state probabilities for the HMM model
!
!
SUBROUTINE get_viterbi_path(K,M,Q,N,n_tau,n_scale,x,z,states,noise,init_t_range,u,phi,t_tau,t_scale,best_path)

USE OMP_LIB
INTEGER, INTENT(IN) :: n_tau,n_scale,M,Q,N,K
INTEGER, INTENT(IN) :: states(M*Q,2)
INTEGER, INTENT(OUT) :: best_path(N,K)
INTEGER :: i, j, kk,ii
INTEGER(KIND=4) :: S
INTEGER :: target_ind, target_m, target_q,init_t_range
INTEGER :: tau_step,scale_step,transition_inds(M*Q*n_tau*3,2),states_0(0:(M*Q),2)
INTEGER :: tar_inds(n_tau*3)
DOUBLE PRECISION, INTENT(IN) :: t_tau(n_tau), t_scale(n_scale)
DOUBLE PRECISION, INTENT(IN) :: x(N,K),z(M),noise,u(K),phi(Q)
DOUBLE PRECISION :: V_mat(M*Q,N,K), possible_paths(M*Q,N,K), log_likelihood(n_tau*3,K)
DOUBLE PRECISION :: transition_prob(M*Q*n_tau*3),prob_sum
DOUBLE PRECISION :: state_p(n_tau*3,K)


!!!!!!!!!!!!!!!!!!!!!!!!!
! Set some constants and initialize arrays to zero
!
! Number of states
S = M*Q
! Number of states that can be transitioned to from a given state
ntars = 3*n_tau
! Initialize log_likelihood to 0 to V_mat to all -Infinity
log_likelihood(:,:) = 0.d0
V_mat(:,:,:) = LOG(log_likelihood(1,1))
! states_0 extends the states array to include the "dump" state for dealing
!     with boundary transitions
states_0(0,:) = [1,1]
states_0(1:S,:) = states

!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the transition matrix of the HMM
!     Implementation attempts to reduce computations 
!     when transition matrix is sparse.
DO i=1,S
    ! Fill in all states that can be reached from state i
    ! NOTE: It is inefficient to check the up and down scale steps
    !           when there is only 1 scale state. 
    DO tau_step=1,n_tau; DO scale_step=-1,1
        ! Obtain the transition probabilities for the scale and time transitions
        p_tau=t_tau(tau_step)
        IF (ABS(scale_step)==1) THEN 
            p_scale=t_scale(2)
        ELSE
            p_scale=t_scale(1)
        ENDIF
        target_m = states(i,1)+tau_step
        ! Check to see if target m is passed the last latent time
        IF (target_m .GT. M) THEN
            p_tau=0.
            ! This is a little sloppy - prevents array out of bounds later 
            !     State 0 is always set to be a zero probability state that 
            !     is the target of all out of state space steps
            target_m=0
        ENDIF
        target_q = states(i,2)+scale_step
        ! Check to see if target q is outside the scale range
        IF (target_q .LT. 1 .OR. target_q .GT. Q) THEN
            p_scale=0.
            ! This is sloppy - See comment above about doing this in tau step
            target_q=0
        ENDIF
        IF (target_m == 0 .OR. target_q == 0) THEN
            target_ind = 0
        ELSE
            target_ind = target_m+(target_q-1)*M
        ENDIF
        ! Set transition probability for jump from i to target_ind 
        !     Also set matrix to carry the inds, first column holds current index,
        !     second column holds target index. 
        ! NOTE: transition_inds in principle only needs to have one column, 
        !     the current index defines the ordering (first ntars values are i, next are i+1
        !     and so on). Easily could be dropped to save memory. I'm leaving it in for now
        !     for clarity of what transition_inds is doing.
        transition_prob((i-1)*ntars+(tau_step-1)*3+scale_step+2) = p_scale*p_tau
        transition_inds((i-1)*ntars+(tau_step-1)*3+scale_step+2,2) = target_ind
        transition_inds((i-1)*ntars+(tau_step-1)*3+scale_step+2,1) = i
    ENDDO; ENDDO
    ! Normalize each set
    prob_sum = SUM(transition_prob(((i-1)*ntars+1):(i*ntars)))
    if(prob_sum .GT. 0) transition_prob(((i-1)*ntars+1):(i*ntars)) = &
                                    transition_prob(((i-1)*ntars+1):(i*ntars))/prob_sum
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!
! Run forward-backward algoritm for calculating the state probabilities
! k iterates over the replicate time series
!
! OpenMP parallelism for shared memory, multi-core parallelization
!     Be careful of race conditions if you edit this!!!
!
!$OMP PARALLEL DO DEFAULT(SHARED), FIRSTPRIVATE(i,j,num_tars,tar_inds,ii)
DO kk=1,K
    ! Initialize the forward and backward probability matrices
    DO i=1,Q
        ! NOTE: NEED TO GET RID OF THIS HARDCODED 100 HERE
        !    Its here to allow the first time of the day to vary more around
        !    the profile than subsequent steps.
        IF (x(1,kk) .LT. -1.d+24) THEN
            V_mat((1+(i-1)*M):(init_t_range+(i-1)*M),1,kk) = LOG(1.0d0/(1.d0*init_t_range*Q))
        ELSE
            V_mat((1+(i-1)*M):(init_t_range+(i-1)*M),1,kk) = LOG(1.0d0/(1.d0*init_t_range*Q)) - &
                    (x(1,kk) - z(1:init_t_range)*u(kk)*phi(i))**2/(2.0d0*noise)
        ENDIF
    ENDDO
    ! Iterate over observations within a series
    DO i=2,N
        ! For each target state, find the current state that maximizes the probability of
        !     reaching that state. The loop over j below iterates over CURRENT states in 
        !     order to avoid a search over the target state vector that slowed the code
        !     significantly. For each current state, check to see if it maximizes the 
        !     probability of reaching it's possible target states (of all the current 
        !     states checked so far). Avoiding the search is really crucial with many states
        DO j=1,S
            tar_inds = transition_inds(((j-1)*ntars+1):(j*ntars),2)
            IF (x(i,kk) .LT. -1.d+24) THEN
                log_likelihood(:,kk) = 0.0
            ELSE
                log_likelihood(:,kk) = &
                        -(x(i,kk) - z(states(tar_inds,1))*u(kk)*phi(states(tar_inds,2)))**2/(2.0d0*noise)
            ENDIF
            state_p(:,kk) = LOG(transition_prob(((j-1)*ntars+1):(j*ntars))) + &
                    V_mat(j,i-1,kk) + &
                    log_likelihood(:,kk)
            DO ii=1,ntars
                IF ( state_p(ii,kk) .GT. V_mat(tar_inds(ii),i,kk) ) THEN
                    V_mat(tar_inds(ii),i,kk) = state_p(ii,kk)
                    possible_paths(tar_inds(ii),i-1,kk) = j
                END IF
            ENDDO
        ENDDO
    ENDDO
    ! Find final state with maximum probability and work backwards to get the full path
    best_path(N,kk) = MAXLOC(V_mat(:,N,kk),1)
    DO i=N-1,1,-1
        best_path(i,kk) = possible_paths(best_path(i+1,kk),i,kk)
    ENDDO
ENDDO
!$OMP END PARALLEL DO

END
