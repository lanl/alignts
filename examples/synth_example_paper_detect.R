library(alignts)
library(lhs)
load(file = "alignment_example.RData")

n_sims         <- 500
lhs_size       <- maximinLHS(n_sims,2)
perturb_params <- data.frame(time   = sample(c(4,28),n_sims,replace = TRUE),
                             height = 6*lhs_size[,1]-3)
perturb_params$duration <- 0
for (i in 1:n_sims) perturb_params$duration[i] <- ceiling(lhs_size[i,2]*16)

for (i in 1:n_sims) {
  test_obs [[i]] <- generate_replicate(rep(1/3,3),1,1)
  perturb        <- c(rep(0,perturb_params$time[i]-1),
                      rep(perturb_params$height[i],perturb_params$duration[i]),
                      rep(0,48 - perturb_params$time[i] - perturb_params$duration[i]+1))
  test_y[[i]]  <- test_obs[[i]][,2] + perturb
}

train_resid  <- list()

for (ind in 1:60){ 
  ttt <- proc.time()
  train_resid[[ind]] <- init_particle_filter(200,48,1:3,1:3)
  for (curr_time in 1:48){
    train_resid[[ind]] <- particle_align(t           = curr_time,
                                         x           = train_y[[ind]][curr_time],
                                         aligned_obj = a_y,
                                         particles   = train_resid[[ind]],
                                         ind         = curr_time+1,
                                         periodic    = TRUE)
  }
  print(proc.time()-ttt)
}

mean_residuals  <- lapply(train_resid,
                          function(x) apply(sapply(x[1:200],function(y)y$residual[-1]),1,weighted.mean,w=x$w))
cuts_train  <- find_variance_partitions(lapply(1:30, function(x)1:48),
                                        mean_residuals,
                                        n_long   = 48,
                                        h        = 4,
                                        crit     = 1.5,
                                        max_cuts = 20)

adj_mean_resids <- mapply(adjust_residuals,
                          x = mean_residuals, 
                          t = lapply(1:30, function(x)1:48), 
                          MoreArgs = list(cuts=cuts_train), SIMPLIFY = FALSE)

cusum_train_mean = mean(unlist(adj_mean_resids))


#####
# Lists to hold online residuals for each "current time"
#    as if the new times have not yet been observed
#
test_resid  <- list()

for (ind in 1:n_sims){ 
  ttt <- proc.time()
  test_resid[[ind]] <- init_particle_filter(200,48,1:3,1:3)
  for (curr_time in 1:48){
    test_resid[[ind]] <- particle_align(t           = curr_time,
                                        x           = test_y[[ind]][curr_time],
                                        aligned_obj = a_y,
                                        particles   = test_resid[[ind]],
                                        ind         = curr_time+1,
                                        periodic    = TRUE)
  }
  print(proc.time()-ttt)
}

detect_stat      <- matrix(0,n_sims,48)
detect_stat_test <- matrix(0,n_sims,48)

detect_stat_clean <- list()

max_stat = 4.1
for (ind in 1:length(adj_mean_resids)){
  
  t   <- 1:48
  
  detect_test = cusum(adj_mean_resids[[ind]],sqrt(a_y$noise),cusum_train_mean)
  detect_stat_clean[[ind]] <- detect_test
}

print(mean(sapply(detect_stat_clean, function(x) sum(x[-1] == 0) )))

plot(0,0,
     type = "n",
     xlim = c(0,48),
     ylim = range(unlist(detect_stat_clean)))
for (ii in 1:30) lines(1:48,
                       detect_stat_clean[[ii]],
                       lwd = 2,
                       col = rgb(1-(ii-1)/30,0,(ii-1)/30,0.6))


mean_residuals_d  <- lapply(test_resid,
                            function(x)apply(sapply(x[1:200],function(y)y$residual[-1]),1,weighted.mean,w=x$w))
adj_mean_resids_d <- mapply(adjust_residuals,
                            x = mean_residuals_d, 
                            t = lapply(1:500, function(x)1:48), 
                            MoreArgs = list(cuts=cuts_train), SIMPLIFY = F)

detect_stat_dirty <- list()

for (ind in 1:length(test_y)){
  
  t   <- 1:48
  y   <- test_y[[ind]]
  
  detect_test = cusum(adj_mean_resids_d[[ind]],sqrt(a_y$noise),cusum_train_mean)
  detect_stat_dirty[[ind]] <- detect_test
}

print(mean(sapply(detect_stat_dirty, function(x) sum(x[-1] == 0) )))

plot(0,0,
     type = "n",
     xlim = c(0,48),
     ylim = range(unlist(detect_stat_dirty)))
for (ii in 1:500) lines(1:48,
                        detect_stat_dirty[[ii]],
                        lwd = 2,
                        col = rgb(1-(ii-1)/500,0,(ii-1)/500,0.6))

pdf("synth_filtered_residuals_dirty.pdf",width=10,height=5)
par(mar = c(5, 5, 1, 1))
plot(0,0,
     type = "n",
     xaxt = "n",
     ylab = "Filtered Residuals",
     xlab = "Observed Time",
     xlim = c(0,48),
     ylim = c(-5,5),
     cex.axis = 1.5,
     cex.lab = 1.8)
axis(1,at=c(1,12,24,36,48),
     labels=c("12 AM","6 AM", "12 PM", "6 PM", "12 AM"),
     cex.axis=1.5)
for (ii in 1:500) lines(1:48,
                        adj_mean_resids_d[[ii]],
                        lwd = 2,
                        col = rgb(1-(ii-1)/500,0,(ii-1)/500,0.6))
dev.off()


pdf("synth_display_detections.pdf",width=10,height=4)
par(mar = c(5, 5, 1, 1))
false_pos = sapply(detect_stat_dirty, function(x) sum(x[-1] == 0) )==0
plot(perturb_params$duration[!false_pos]*1440/48,
     perturb_params$height[!false_pos],
     xlab = "Anomaly Duration (minutes)",
     ylab = "Anomaly Magnitude (C)",
     cex = 2.0,
     cex.lab = 1.7,
     cex.axis = 1.5,
     ylim = c(-3.2,3.2),
     yaxt="n",
     pch=19,
     col="dark blue")
points(perturb_params$duration[false_pos]*1440/48,
       perturb_params$height[false_pos],
       xlab = "Anomaly Duration (minutes)",
       ylab = "Anomaly Magnitude (C)",
       cex = 1.5,
       cex.lab = 1.7,
       cex.axis = 1.5,
       pch=19,
       col="grey")
axis(2,at=c(-3:3),
     labels=c("-3","","","0","","","3"),
     cex.axis=1.5)
abline(h=0,lty=2, lwd=4, col=rgb(0,0,0,0.4))
dev.off()

plot(0,0,
     type = "n",
     xlim = c(0,48),
     ylim = range(unlist(test_y)))
for (ii in 1:500) lines(1:48,
                        test_y[[ii]],
                        lwd = 2,
                        col = rgb(1-(ii-1)/500,0,(ii-1)/500,0.6))



test_obs_clean <- list()
test_y_clean   <- list()

n_sims <- 500

for (i in 1:n_sims) {
  test_obs_clean[[i]] <- generate_replicate(rep(1/3,3),1,1)
  test_y_clean[[i]]  <- test_obs[[i]][,2]
}

test_resid_clean  <- list()

for (ind in 1:n_sims){ 
  test_resid_clean[[ind]] <- init_particle_filter(200,48,1:3,1:3)
  for (curr_time in 1:48){
    test_resid_clean[[ind]] <- particle_align(t           = curr_time,
                                              x           = test_y_clean[[ind]][curr_time],
                                              aligned_obj = a_y,
                                              particles   = test_resid_clean[[ind]],
                                              ind         = curr_time+1,
                                              periodic    = TRUE)
  }
}

mean_residuals_clean  <- lapply(test_resid_clean,
                                function(x)apply(sapply(x[1:200],function(y)y$residual[-1]),1,weighted.mean,w=x$w))
adj_mean_resids_clean <- mapply(adjust_residuals,
                                x = mean_residuals_clean, 
                                t = lapply(1:500, function(x)1:48), 
                                MoreArgs = list(cuts=cuts_train), SIMPLIFY = F)

detect_new_clean <- list()

for (ind in 1:length(test_y_clean)){
  
  t   <- 1:48
  y   <- test_y_clean[[ind]]
  
  detect_test = cusum(adj_mean_resids_clean[[ind]],sqrt(a_y$noise),cusum_train_mean)
  detect_new_clean[[ind]] <- detect_test
}

print(mean(sapply(detect_new_clean, function(x) sum(x[-1] == 0) )))

pdf("synth_filtered_residuals_clean.pdf",width=10,height=5)
par(mar = c(5, 5, 1, 1))
plot(0,0,
     type = "n",
     xaxt = "n",
     ylab = "Filtered Residuals",
     xlab = "Observed Time",
     xlim = c(0,48),
     ylim = c(-5,5),
     cex.axis = 1.5,
     cex.lab = 1.8)
axis(1,at=c(1,12,24,36,48),
     labels=c("12 AM","6 AM", "12 PM", "6 PM", "12 AM"),
     cex.axis=1.5)
for (ii in 1:500) lines(1:48,
                        adj_mean_resids_clean[[ii]],
                        lwd = 2,
                        col = rgb(1-(ii-1)/500,0,(ii-1)/500,0.6))
dev.off()


pdf("synth_filtered_residuals_both.pdf",width=10,height=6)
#par(mfrow=c(2,1))
layout(matrix(c(1,2),2,1), widths=1, heights= c(1,1.35))
par(mar = c(0, 5, 1, 1))
plot(0,0,
     type = "n",
     xaxt = "n",
     ylab = "Residuals",
     xlab = "",
     xlim = c(0,48),
     ylim = c(-5,5),
     xaxt = "n",
     yaxt = "n",
     cex.axis = 1.5,
     cex.lab = 1.7)
axis(2,at=c(-4,0,4),
     labels=c("-4","0", "4"),
     cex.axis=1.5)
for (ii in 1:500) lines(1:48,
                        adj_mean_resids_clean[[ii]],
                        lwd = 2,
                        col = rgb(1-(ii-1)/500,0,(ii-1)/500,0.6))

par(mar = c(5, 5, 0, 1))
plot(0,0,
     type = "n",
     xaxt = "n",
     yaxt = "n",
     ylab = "Residuals",
     xlab = "Observed Time",
     xlim = c(0,48),
     ylim = c(-5,5),
     cex.axis = 1.5,
     cex.lab = 1.7)
axis(2,at=c(-4,0,4),
     labels=c("-4","0", "4"),
     cex.axis=1.5)
axis(1,at=c(1,12,24,36,48),
     labels=c("12 AM","6 AM", "12 PM", "6 PM", "12 AM"),
     cex.axis=1.5)
for (ii in 1:500) lines(1:48,
                        adj_mean_resids_d[[ii]],
                        lwd = 2,
                        col = rgb(1-(ii-1)/500,0,(ii-1)/500,0.6))
dev.off()

