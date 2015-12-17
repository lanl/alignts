!------------------------------------------------
! Function for calculating the state probabilities for the HMM model
!
!
SUBROUTINE e_step_periodic(K,M,Q,N,n_tau,n_scale,x,z,states,noise,init_t_range,u,phi,t_tau,t_scale,state_prob)

USE OMP_LIB
INTEGER, INTENT(IN) :: n_tau,n_scale,M,Q,N,K
INTEGER :: i, j, kk, ii,iii
INTEGER(KIND=4) :: S
INTEGER :: target_ind, target_m, target_q, init_t_range
INTEGER, INTENT(IN) :: states(M*Q,2)
INTEGER :: tau_step,scale_step,ntars,transition_inds(M*Q*n_tau*3,2),states_0(0:(M*Q),2)
DOUBLE PRECISION, INTENT(IN) :: t_tau(n_tau), t_scale(n_scale)
DOUBLE PRECISION, INTENT(IN) :: x(N,K),z(M),noise,u(K),phi(Q)
DOUBLE PRECISION, INTENT(OUT) :: state_prob(M*Q,N,K)
DOUBLE PRECISION, PARAMETER :: pi = ACOS(-1.d0)
DOUBLE PRECISION :: forward_mat(N+1,0:(M*Q),K), backward_mat(0:(M*Q),N+1,K), obs_likelihood(n_tau*3,K)
DOUBLE PRECISION :: transition_prob(M*Q*n_tau*3),prob_sum, normC, p_scale, p_tau

!!!!!!!!!!!!!!!!!!!!!!!!!
! Set some constants and initialize arrays to zero
!
! Set Gaussian normalizing constant
normC = 1.0d0/SQRT(2.0d0*pi*noise)
! Number of states
S = M*Q
! Number of states that can be transitioned to from a given state
ntars = 3*n_tau
! Initialize arrays to zero
transition_prob(:) = 0.d0
obs_likelihood(:,:) = 0.0d0
forward_mat(:,:,:) = 0.0d0
backward_mat(:,:,:) = 0.0d0
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
            ! Target beginning state with periodicity
            target_m = target_m - M
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
! k iterates over the time series
!
! OpenMP parallelism for shared memory, multi-core parallelization
!
!$OMP PARALLEL DO DEFAULT(SHARED), FIRSTPRIVATE(i,j,ii,iii,prob_sum)
DO kk=1,K
    ! Initialize the forward and backward probability matrices
    DO i=1,Q
        forward_mat(1,(1+(i-1)*M):(init_t_range+(i-1)*M),kk) = 1.0d0/(1.d0*init_t_range*Q)
    ENDDO
    prob_sum = SUM(forward_mat(1,1:S,kk))
    if(prob_sum .GT. 0) forward_mat(1,1:S,kk) = forward_mat(1,1:S,kk)/prob_sum
    backward_mat(:,N+1,kk) = 1.0d0
    ! Iterate over observations in each replicate series
    DO i=1,N
        ! The value here was a workaround because FORTRAN didn't like being passed
        !     NA values from R, so they were replaced with hugely tiny values. For 
        !     times when there was no observation, "emission" matrix can be treated
        !     as an identity, so forward probabilities are just the current times the 
        !     transition probabilities.
        IF (x(i,kk) .LT. -1.d+24) THEN
            DO j=1,S
                ii  = (j-1)*ntars+1
                iii = (j-1)*ntars+ntars
                forward_mat(i+1,transition_inds(ii:iii,2),kk) = &
                        forward_mat(i+1,transition_inds(ii:iii,2),kk) + &
                        forward_mat(i,j,kk)*transition_prob(ii:iii)
            ENDDO
            prob_sum = SUM(forward_mat(i+1,1:S,kk))
            if(prob_sum .GT. 0) forward_mat(i+1,1:S,kk) = forward_mat(i+1,1:S,kk)/prob_sum
        ELSE 
            DO j=1,S
                ii  = (j-1)*ntars+1
                iii = (j-1)*ntars+ntars
                obs_likelihood(:,kk) = &
                        normC*EXP(-(x(i,kk) - z(states(transition_inds(ii:iii,2),1))* & 
                        u(kk)*phi(states(transition_inds(ii:iii,2),2)))**2/(2.0d0*noise))
                forward_mat(i+1,transition_inds(ii:iii,2),kk) = &
                        forward_mat(i+1,transition_inds(ii:iii,2),kk) + &
                        forward_mat(i,j,kk)*obs_likelihood(:,kk)*transition_prob(ii:iii)
            ENDDO
            prob_sum = SUM(forward_mat(i+1,1:S,kk))
            if(prob_sum .GT. 0) forward_mat(i+1,1:S,kk) = forward_mat(i+1,1:S,kk)/prob_sum
        ENDIF
        IF (x(N-i+1,kk) .LT. -1.d+24) THEN
            DO j=1,S
                ii  = (j-1)*ntars+1
                iii = (j-1)*ntars+ntars
                backward_mat(j,N-i+1,kk) = DOT_PRODUCT(transition_prob(ii:iii), & 
                        backward_mat(transition_inds(ii:iii,2),N-i+2,kk))
            ENDDO
            prob_sum = SUM(backward_mat(1:S,N-i+1,kk))
            if(prob_sum .GT. 0) backward_mat(1:S,N-i+1,kk) = backward_mat(1:S,N-i+1,kk)/prob_sum
        ELSE
            ! For each state, calculate backward probability
            DO j=1,S
                ii  = (j-1)*ntars+1
                iii = (j-1)*ntars+ntars
                obs_likelihood(:,kk) = &
                        normC*EXP(-(x(N-i+1,kk) - z(states(transition_inds(ii:iii,2),1))* &
                        u(kk)*phi(states(transition_inds(ii:iii,2),2)))**2/(2.0d0*noise))
                backward_mat(j,N-i+1,kk) = & 
                        DOT_PRODUCT(transition_prob(ii:iii), & 
                                obs_likelihood(:,kk)*backward_mat(transition_inds(ii:iii,2),N-i+2,kk))
            ENDDO
            prob_sum = SUM(backward_mat(1:S,N-i+1,kk))
            if(prob_sum .GT. 0) backward_mat(1:S,N-i+1,kk) = backward_mat(1:S,N-i+1,kk)/prob_sum
        ENDIF
    ENDDO
    ! Calculate state probabilities for all times from forward and backward probabilities
    state_prob(:,:,kk) = TRANSPOSE(forward_mat(2:(N+1),1:S,kk))*backward_mat(1:S,2:(N+1),kk)
    DO i=1,N
        prob_sum = SUM(state_prob(:,i,kk))
        if(prob_sum .GT. 0) state_prob(:,i,kk) = state_prob(:,i,kk)/prob_sum        
    ENDDO
ENDDO
!$OMP END PARALLEL DO

END
