library(lhs)
library(alignts)
source("source_functions.R")

####
# Generate training and validation data
#
train_obs <- list()
train_y   <- list()

for (i in 1:60) {
  train_obs [[i]] <- generate_replicate(rep(1/3,3),1,1)
  train_y[[i]]    <- train_obs[[i]][,2]
}

valid_obs <- list()
valid_y   <- list()

for (i in 1:30) {
  valid_obs [[i]] <- generate_replicate(rep(1/3,3),1,1)
  valid_y[[i]]    <- valid_obs[[i]][,2]
}


#####
# Generate test data, some with perturbations
#
test_obs <- list()
test_y   <- list()

#####
# Fit alignment model with training data
#
bench_sse <- rep(Inf,10)
lambdas   <- 10^seq(2,-2,length=10)
ta_y      <- list()
tt        <- 1:48

for (ind in 1:10){
  ta_y[[ind]] <- align_series_EM(t                = tt,
                                 x                = train_y,
                                 J                = 3,
                                 periodic         = TRUE,
                                 lambda           = lambdas[ind],
                                 iters            = 400,
                                 init_z           = apply( sapply(train_y,function(x)x), 1, mean),
                                 parallel_cores   = 1,
                                 scales           = c(1,1.1,1.2), 
                                 transition_scale = c(0.8,0.2),
                                 transition_tau   = dbinom(0:(3-1),3-1,0.5)
  )
  p_y <- predict_series_EM(t                = tt,
                           x                = valid_y,
                           aligned_obj      = ta_y[[ind]],
                           J                = 3,
                           periodic         = TRUE,
                           transition_scale = c(0.8,0.2),
                           transition_tau   = dbinom(0:(3-1),3-1,0.5),
                           iters            = 400,
                           parallel_cores   = 1
  )
  bench_sse[ind] <- sum(unlist(get_residuals(t = tt, x = valid_y, aligned_obj = p_y))^2)
}

ind    <- which(bench_sse < 1.05*min(bench_sse))[1]
lambda <- lambdas[ind]
print(bench_sse)

a_y <- align_series_EM(t                = tt,
                       x                = train_y,
                       J                = 3,
                       periodic         = TRUE,
                       lambda           = lambda,
                       iters            = 400,
                       init_z           = apply( sapply(train_y,function(x)x), 1, mean),
                       parallel_cores   = 1,
                       scales           = c(1,1.1,1.2), 
                       transition_scale = c(0.8,0.2),
                       transition_tau   = dbinom(0:(3-1),3-1,0.5)
)

save(a_y,tt,train_y, valid_y, test_y, train_obs, valid_obs, test_obs,file = "alignment_example.RData")

pdf("synth_estimated_curve.pdf",width=10,height=8)
x_scaled <- train_y
for (i in 1:n_series) x_scaled[[i]] <- train_y[[i]]/a_y$u[i]/a_y$scales[states$q[a_y$viterbi_path[t,i]]]
par(mar = c(6, 6, 2, 2))
t <- 1:48
n_series <- 60
states <- a_y$states
tau <- a_y$viterbi_tau
plot(rep(1:48, n_series), unlist(x_scaled), type = "n", xaxs="i",yaxs="i",xlim=c(0,50),
     ylab = "Scaled, Aligned Output", xlab = "Time", main = "", cex.lab = 1.5, 
     cex.axis = 1.5,ylim=c(2,8),xaxt="n")
axis(1,at=c(1,12,24,36,48),labels=c("12 AM","6 AM", "12 PM", "6 PM", "12 AM"),cex.axis=1.5)
cols <- rgb((1:n_series - 1)/n_series, 0, 1 - (1:n_series - 
                                                 1)/n_series, 0.25)
for (i in 1:n_series) points(tau[t, i], x_scaled[[i]], 
                             pch = 19, cex = 0.7, col = cols[i], type = "p", 
                             lwd = 3)
lines(a_y$tau, a_y$z, type = "l", lwd = 3,col="green")
lines(1:96/2,true_curve(1:96),lwd=3)
dev.off()

##### Remake Residual plot
par(mar = c(5, 5, 2, 2))
residuals <- get_residuals(1:48,train_y,a_y)
n_series <- length(residuals)
cols <- rgb((seq_along(residuals) - 1)/n_series, 0, 1 - (seq_along(residuals) -
                                                           1)/n_series, 0.6)
plot(rep(t, n_series), unlist(residuals), type = "n",
     ylab = "Model Residuals", xlab = "Time", cex.axis = 1.7, cex.lab = 1.7,xaxt="n")
axis(1,at=c(1,12,24,36,48),labels=c("12 AM","6 AM", "12 PM", "6 PM", "12 AM"),cex.axis=1.7)
dummy <- lapply(1:n_series, function(i, r, t, cols) {
  lines(t, r[[i]], col = cols[i], lwd = 2)
}, r = residuals, t = t, cols = cols)


par(mfrow=c(1,2))
par(mar=c(5,5,3,2))
qqnorm(unlist(residuals),pch=19,cex=0.75,cex.lab=1.7,cex.axis=1.7)
qqline(unlist(residuals),lwd=2)
acf10 <- sapply(0:10,acf_reps, x=residuals)
plot(0:10,acf10,type="h", lwd=4, xlab="Lag",ylab="ACF",
     cex.lab=1.7,cex.axis=1.7)
abline(h=0)

par(mfrow=c(1,1))
par(mar = c(4.5, 4.5, 3, 3))

plot_aligned_x(tt,train_y,a_y)
