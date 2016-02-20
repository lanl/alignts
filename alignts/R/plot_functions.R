##########
# Function definitions for various plotting of the alignment output
#
# Required libraries:
#   NONE
#
##########

#####
# Plot the raw (x,t) data from the unevenly sampled lists
#
plot_raw <- function(t,x,lwd=2,xlab="",ylab=""){
  n_series <- length(x)
  if(is.list(t)){
    plot(unlist(t),unlist(x),type="n",xlab=xlab,ylab=ylab,cex.lab=1.5,cex.axis=1.5)
    for (i in 1:n_series) {
      lines(t[[i]],x[[i]],lwd=lwd,col=rgb((i-1)/n_series,0,1-(i-1)/n_series,0.6))
    }
  } else{
    plot(rep(t,n_series),unlist(x),type="n",xlab=xlab,ylab=ylab,cex.lab=1.5,cex.axis=1.5)
    for (i in 1:n_series) {
      lines(t,x[[i]],lwd=lwd,col=rgb((i-1)/n_series,0,1-(i-1)/n_series,0.6))
    }
  }
}

#####
# Plot the scaled, time-warped observed output curve for each replicate 
#     series. This shows the approximate alignment with the variability 
#     around the collective latent profile being the scaled noise.
#
plot_aligned_x <- function(t,x,aligned_obj,type="p",xlab="Latent Time",ylab="Scaled Observed \n Output",main="Viterbi Path"){
  x_scaled <- x
  n_series <- length(x)
  states <- aligned_obj$states
  tau <- aligned_obj$viterbi_tau
  if (is.list(t)){
    for (i in 1:n_series) x_scaled[[i]]<- x[[i]]/aligned_obj$u[i]/aligned_obj$scales[states$q[aligned_obj$viterbi_path[t[[i]],i]]]
    par(mar=c(6,6,2,2))
    plot(unlist(t),unlist(x_scaled),type="n",ylab=ylab,xlab=xlab,
         main=main,cex.lab=1.4,cex.axis=1.4,cex.main=1.35)
    cols <- rgb((1:n_series-1)/n_series,0,1-(1:n_series-1)/n_series,0.4)
    for (i in 1:n_series) points(tau[t[[i]],i],x_scaled[[i]],pch=19,cex=0.7,col=cols[i],type=type,lwd=3)
    lines(aligned_obj$tau,aligned_obj$z,type="l",lwd=3)
    par(mar=c(4.5,4.5,3,3))
  }else{
    for (i in 1:n_series) x_scaled[[i]]<- x[[i]]/aligned_obj$u[i]/aligned_obj$scales[states$q[aligned_obj$viterbi_path[t,i]]]
    par(mar=c(6,6,2,2))
    plot(rep(t,n_series),unlist(x_scaled),type="n",ylab=ylab,xlab=xlab,
         main=main,cex.lab=1.4,cex.axis=1.4,cex.main=1.35)
    cols <- rgb((1:n_series-1)/n_series,0,1-(1:n_series-1)/n_series,0.4)
    for (i in 1:n_series) points(tau[t,i],x_scaled[[i]],pch=19,cex=0.7,col=cols[i],type=type,lwd=3)
    lines(aligned_obj$tau,aligned_obj$z,type="l",lwd=3)
    par(mar=c(4.5,4.5,3,3))
  }
}

#####
# Plot an observed time series with it's alignment to the characteristic
#     latent profile. Useful for illustrating exactly what the method is doing
#
warp_example <- function(t, x, aligned_obj, series_ind, plot_every=1){
  # v_inds is the vector of viterbi states for the HMM
  n <- length(x[[series_ind]])
  plot_inds <- seq(1,n,by=plot_every)
  x <- x[[series_ind]][plot_inds]
  if(is.list(t)){
    t <- t[[series_ind]][plot_inds]
  } else{
    t <- t[plot_inds]    
  }
  z <- aligned_obj$z
  tau <- aligned_obj$tau
  M <- length(z)
  v_inds <- (aligned_obj$viterbi_path[,series_ind]-1) %% M +1
  plot(t,x,pch=19,cex=0.8,xlab="Time",ylab="Output", main=paste("Series #",series_ind,sep=""),
       xlim=range(tau),ylim=range(c(x,z)),col=rgb(1,0,0,0.75),cex.lab=1.4,cex.axis=1.4)
  lines(tau,z,lwd=3)
  for(i in 1:length(x)) segments(x0 = t[i],y0=x[i],
                                 x1=tau[v_inds[t[i]]],y1=z[v_inds[t[i]]],
                                 lwd=1,col=rgb(1,0,0,0.5))
}

#####
# Plot an observed time series with it's alignment to the characteristic
#     latent profile, but in steps so that each parameter's effect can be
#     views: global scaling, then local scaling, then time warping.
#
warp_steps_illustrate <- function(t, x, aligned_obj, series_ind, plot_every=1, plot_separate=FALSE){
  n <- length(x[[series_ind]])
  plot_inds <- seq(1,n,by=plot_every)
  z <- aligned_obj$z
  tau <- aligned_obj$tau
  u <- aligned_obj$u[series_ind]
  states <- aligned_obj$states
  M <- length(z)
  v_inds <- (aligned_obj$viterbi_path[plot_inds,series_ind]-1) %% M +1
  phi <- aligned_obj$scales[states$q[aligned_obj$viterbi_path[plot_inds,series_ind]]]
  x <- x[[series_ind]][plot_inds]
  if(is.list(t)){
    t <- t[[series_ind]][plot_inds]
  } else{
    t <- t[plot_inds]    
  }
  point_size = 1.0
  text_size = 1.4
  line_size = 2
  curve_size = 4
  if(plot_separate){
    par(mfrow=c(1,1))
  } else{
    par(mfrow=c(2,2))
  }
  # Start with the full fit, leaving only residual variability
  #
  plot(t,x,pch=19,cex=point_size,xlab="Time",ylab="Output", cex.lab=text_size, cex.axis=text_size,
       xlim=range(tau),ylim=range(c(x,z)),col=rgb(1,0,0,0.75), main="Model Fit Showing Only Residual Noise")
  lines(tau,z,lwd=curve_size)
  points(t,z[v_inds]*u*phi,pch=19,cex=point_size,col=rgb(0.5,0,1,0.75))
  for(i in 1:length(x)) segments(x0 = t[i],y0=x[i],
                                 x1=t[i],y1=z[v_inds[i]]*u*phi[i],
                                 lwd=line_size,col=rgb(1,0,0,0.5))

  # Remove time-warping to put the hidden state values back on latent time
  #
  plot(t,x,pch=19,cex=point_size,xlab="Time",ylab="Output", cex.lab=text_size, cex.axis=text_size,
       xlim=range(tau),ylim=range(c(x,z)),col=rgb(1,0,0,0.75), main="Remove Time-Warping")
  lines(tau,z,lwd=curve_size)
  points(tau[v_inds],z[v_inds]*u*phi,pch=19,cex=point_size,col=rgb(0.5,0,1,0.75))
  for(i in 1:length(x)) segments(x0 = t[i],y0=x[i],
                                 x1=tau[v_inds[i]],y1=z[v_inds[i]]*u*phi[i],
                                 lwd=line_size,col=rgb(1,0,0,0.5))

  # Remove local scaling, leaving only uniform scaling
  #
  plot(t,x,pch=19,cex=point_size,xlab="Time",ylab="Output", cex.lab=text_size, cex.axis=text_size,
       xlim=range(tau),ylim=range(c(x,z)),col=rgb(1,0,0,0.75), main="Remove Local Scaling")
  lines(tau,z,lwd=curve_size)
  points(tau[v_inds],z[v_inds]*u,pch=19,cex=point_size,col=rgb(0.5,0,1,0.75))
  for(i in 1:length(x)) segments(x0 = t[i],y0=x[i],
                                 x1=tau[v_inds[i]],y1=z[v_inds[i]]*u,
                                 lwd=line_size,col=rgb(1,0,0,0.5))

  # Plot observed series aligned to latent series by removing uniform scaling
  #
  plot(t,x,pch=19,cex=point_size,xlab="Time",ylab="Output", cex.lab=text_size, cex.axis=text_size,
       xlim=range(tau),ylim=range(c(x,z)),col=rgb(1,0,0,0.75), main="Remove Global Scaling")
  lines(tau,z,lwd=curve_size)
  points(tau[v_inds],z[v_inds],pch=19,cex=point_size,col=rgb(0.5,0,1,0.75))
  for(i in 1:length(x)) segments(x0 = t[i],y0=x[i],
                                 x1=tau[v_inds[i]],y1=z[v_inds[i]],
                                 lwd=line_size,col=rgb(1,0,0,0.5))
      
}

plot_residuals <- function(t,residuals,xlab="Observed Time",ylab="Model Residuals"){
  n_series <- length(residuals)
  cols <- rgb((seq_along(residuals)-1)/n_series,0,1-(seq_along(residuals)-1)/n_series,0.6)
  if(is.list(t)){
    plot(unlist(t),unlist(residuals),type="n", ylab=ylab,xlab=xlab,
         cex.axis=1.4,cex.lab=1.4)
    mapply(lines,t,residuals,col=cols,lwd=2)
  } else{
    plot(rep(t,n_series),unlist(residuals),type="n", ylab=ylab,xlab=xlab,
         cex.axis=1.4,cex.lab=1.4)
    dummy <- lapply(1:n_series,function(i,r,t,cols){lines(t,r[[i]],col=cols[i],lwd=2)},r=residuals,t=t,cols=cols)
  }
}
