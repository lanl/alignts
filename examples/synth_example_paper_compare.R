library(Rlof)

clean_data_y = matrix(unlist(test_y_clean),nrow=48)
dirty_data_y = matrix(unlist(test_y),nrow=48)

all_test = rbind(t(clean_data_y),t(dirty_data_y))

mean_resid_clean = sapply(test_resid_clean,function(x) apply(sapply(x[1:200], function(y) y$residual),1,mean))
mean_resid_dirty = sapply(test_resid,function(x) apply(sapply(x[1:200], function(y) y$residual),1,mean))
mean_resid_clean = mean_resid_clean[2:49,]
mean_resid_dirty = mean_resid_dirty[2:49,]

all_resid = rbind(t(mean_resid_clean),t(mean_resid_dirty))

all_scores = rep(0,1000)
for (ii in 1:1000){
  all_scores[ii] = lof(rbind(all_test[-ii,][1:499,],all_test[ii,]),k=1)[500]
  if (ii %% 20 == 0){
    print(ii)
  }
}

all_scores_resid = rep(0,1000)
for (ii in 1:1000){
  all_scores_resid[ii] = lof(rbind(all_resid[-ii,][1:499,],all_resid[ii,]),k=1)[500]
  if (ii %% 20 == 0){
    print(ii)
  }
}

print(paste0("Mean: ", mean(order(all_scores,decreasing = T)[1:500] > 500)) )
print(paste0("Mean: ", mean(order(all_scores_resid,decreasing = T)[1:500] > 500)) )


centered_scores = rep(0,1000)
for (ii in 1:1000){
  centered_scores[ii] = lof(rbind(smooth_resid[-ii,][1:499,],smooth_resid[ii,]),k=1)[500]
  if (ii %% 20 == 0){
    print(ii)
  }
}
print(paste0("Mean: ", mean(order(centered_scores,decreasing = T)[1:500] > 500)) )

# Compare with Spline Smoothed:

train_data_y = t(matrix(unlist(train_y),nrow=48))


smooth_ss      = smooth.spline(rep(1:48,60), t(train_data_y)[,])
ss_train_resid = train_data_y - t(matrix(predict(smooth_ss,rep(1:48,60))$y,nrow=48))

smooth_resid = all_test - t(matrix(predict(smooth_ss,rep(1:48,1000))$y,nrow=48))


max_stat = 10
detect_smooth = rep(0,60)
for (ind in 1:60){
  t   <- 1:48
  y   <- ss_train_resid[ind, ]
  
  detect_smooth[ind] <- sum( (cusum(y,sd(ss_train_resid),0) == 0)[2:48] )
}
print(mean(detect_smooth))



detect_smooth = rep(0,1000)
for (ind in 1:1000){
  t   <- 1:48
  y   <- smooth_resid[ind, ]
  
  detect_smooth[ind] <- sum( (cusum(y,sd(ss_train_resid),0) == 0)[2:48] )
}
print(mean(detect_smooth[1:500]))

print(mean(detect_smooth[501:1000]))
print(mean(detect_smooth[501:1000] > 0))

