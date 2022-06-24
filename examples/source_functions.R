#
# True latent function
#
#
true_curve <- function(t){
  exp(-1/400*(t-60)^2) * sin(2*pi*t/8+3) + sin(2*pi*t/48) + 4
}

#
# Function to generate replicates given the true curve and model parameters
#
generate_replicate <- function(p_time_step,p_scale_step, global_scaling){
  p_time_step    <- p_time_step/sum(p_time_step)
  global_scaling <- runif(1,0.97,1.03)
  
  times  <- 1
  scales <- c(1.0,1.1,1.2)
  scale  <- sample(1:3,1)
  obs    <- true_curve(1)*global_scaling
  
  for (i in 1:48){
    new_time <- times[length(times)] + sample(1:length(p_time_step),1,prob = p_time_step)
    new_time <- (new_time-1) %% 96 + 1
    scale    <- ifelse(scale %in% c(1,3), sample(c(scale,2),1,prob = c(0.95,0.025)),
                       sample(1:3,1,prob = c(0.025,0.95,0.025)))
    
    times <- c(times,new_time)
    obs   <- c(obs,true_curve(new_time)*global_scaling*scales[scale]+0.3*rnorm(1))
  }
  cbind(1:48*2,obs[-1])
}


#
# Assorted Functions 
# 
acf_reps <- function(x,lag=1){
  if(is.list(x)){
    x         <- lapply(x,function(xx) xx - mean(unlist(x)))
    acf_pairs <- lapply(x,function(xx) cbind(xx[1:(length(xx)-lag)],xx[(1+lag):length(xx)]))
    acf_pairs <- do.call(rbind,acf_pairs)
    denom     <- sum(unlist(x)^2)
  } else{
    x         <- x - mean(x)
    acf_pairs <- cbind(x[1:(length(x)-lag)],x[(1+lag):length(x)])
    denom     <- sum(x^2)
  }
  sum(acf_pairs[,1]*acf_pairs[,2])/denom
}

find_variance_partitions <- function(t, x, n_long = NULL, h = 4, crit = 3.5,max_cuts = 20,make_plots=FALSE){
  n_list <- length(x)
  n_tot  <- length(unlist(x))
  mean_x <- mean(unlist(x))
  if(make_plots) {
    plot(1:1440,1:1440,ylim=range(unlist(x)),type="n")
    lapply(x,function(y) invisible(lines(y,type="l")))
  }
  
  if (is.null(n_long)) n_long <- max(sapply(x,length))
  
  cut_points <- c()
  cut_ratios <- c()
  
  for (i in 1:max_cuts){
    step_r <- rep(0,n_long-2*h+1)
    for (ind in h:(n_long-h) ) {
      before <- unlist(mapply(function(x,t) x[t <  ind], t=t,x=x, SIMPLIFY = F))
      after  <- unlist(mapply(function(x,t) x[t >= ind], t=t,x=x, SIMPLIFY = F))
      
      step_r[ind-h+1] <- sum((after-mean_x)^2)/sum((before-mean_x)^2)*length(before)/length(after)
    }
    
    min_r <- min(step_r)
    max_r <- max(step_r)

    which_min_r <- which.min(step_r) + h
    which_max_r <- which.max(step_r) + h
    
    test_r   <- max(c(max_r, min_r^-1))
    test_ind <- (c(which_max_r,which_min_r))[which.max(c(max_r, min_r^-1))]
    best_r   <- (c(max_r,min_r))[which.max(c(max_r, min_r^-1))]
    
    mean_x   <- mean(unlist(x))
    if (test_r > crit){
      cut_points <- c(cut_points,test_ind)
      cut_ratios <- c(cut_ratios,best_r)
      x <- mapply(function(t,y) c(y[t < test_ind],mean_x+(y[t >= test_ind] - mean_x)/sqrt(best_r)),y=x,t=t, SIMPLIFY = F)
      if(make_plots) {
        plot(1:1440,1:1440,ylim=range(unlist(x)),type="n")
        lapply(x,function(y) invisible(lines(y,type="l")))
        abline(v=test_ind,lwd=2,lty=2,col="red")
      }
      next
    }
    break
  }
  
  list(cut_points = cut_points, cut_ratios = cut_ratios, mean_x = mean_x)
}

adjust_residuals <- function(t, x, cuts){
  for ( i in 1:length(cuts$cut_points) ){
    cut_ind <- min(c(max(t),cuts$cut_points[i]))
    x <- c(x[t < cut_ind],cuts$mean_x+(x[t >= cut_ind] - cuts$mean_x)/sqrt(cuts$cut_ratios[i]))
  }
  
  x
}

cusum <- function(residuals, sd_fit, mean_fit){
  n <- length(residuals)
  
  scaled_resid = (residuals - mean_fit - cumsum(residuals)/1:n)
  scaled_resid = (scaled_resid - cumsum(scaled_resid)/1:n)
  
  c_minus = rep(0,n)
  c_plus  = rep(0,n)
  c_stat  = rep(0,n)
  
  c_minus[1] = min( c(0, scaled_resid[1]) )
  c_plus[1]  = max( c(0, scaled_resid[1]) )
  c_stat[1]  = max( c(abs(c_minus[1]), c_plus[1]) )
  
  for ( ii in 2:n ) {
    c_minus[ii] = min( c(c_minus[ii-1] + scaled_resid[ii], 0) )
    c_plus[ii]  = max( c(c_plus[ii-1]  + scaled_resid[ii], 0) )
    if ( (c_minus[ii] < -max_stat) || (c_plus[ii] > max_stat) ){
      c_minus[ii] = 0
      c_plus[ii]  = 0
    }
    c_stat[ii]  = max( c(abs(c_minus[ii]), c_plus[ii]) )
  }
  
  return(c_stat)
}

