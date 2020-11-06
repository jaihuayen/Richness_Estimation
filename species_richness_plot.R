S <- 200
r <- 500
thres <- 1.5
model <- 'broken'
set.seed(1111)

if(model == 'unif'){
  pi_i <- runif(S)
}else if(model == 'zipf'){
  pi_i <- 1/((1:S)-0.1)
}else if(model == 'broken'){
  pi_i <- rexp(S)
}else if(model == 'homo'){
  pi_i <- rep(1/S,S)
}
pi_norm <- pi_i/sum(pi_i)
cv <- (sd(pi_norm) / mean(pi_norm))

simu_func_Chao <- function(n){
  s_chao <- 0
  for(i in c(1:r)){
    data <- rmultinom(1, n, pi_norm)
    sample_data <- data[data > 0]
    s_obs <- sum(sample_data > 0)
    f1 <- sum(sample_data == 1)
    f2 <- sum(sample_data == 2)
    if(f2 != 0){
      s_chao[i] <- s_obs + (n-1)/ n * (f1^2)/(2*f2)
    }else{
      s_chao[i] <- s_obs + (n-1)/ n * (f1*(f1-1))/2
    }
  }
  return(mean(s_chao))
}
simu_func_proposed <- function(n){
  s_est <- 0
  for(i in c(1:r)){
    data <- rmultinom(1, n, pi_norm)
    sample_data <- data[data > 0]
    p_i <- sample_data / sum(sample_data)
    f1 <- sum(sample_data == 1)
    f2 <- sum(sample_data == 2)
    C <- 1-f1/n
    s_obs <- sum(sample_data > 0)
    if(f2 != 0){
      if(C < thres){
        s_est[i] <- sum(1/(1-(1-p_i)^n)) + (n-1)/ n * (f1^2)/(2*f2)
      }else{
        s_est[i] <- sum(1/(1-(1-p_i)^n))
      }
    }else{
      if(C < thres){
        s_est[i] <- sum(1/(1-(1-p_i)^n)) + (n-1)/ n * (f1*(f1-1))/2
      }else{
        s_est[i] <- sum(1/(1-(1-p_i)^n))
      }
    }
  }
  return(mean(s_est))
}

simu_func_Chao_RMSE <- function(n){
  s_chao <- 0
  for(i in c(1:500)){
    data <- rmultinom(1, n, pi_norm)
    sample_data <- data[data > 0]
    s_obs <- sum(sample_data > 0)
    f1 <- sum(sample_data == 1)
    f2 <- sum(sample_data == 2)
    if(f2 != 0){
      s_chao[i] <- s_obs + (n-1)/ n * (f1^2)/(2*f2)
    }else{
      s_chao[i] <- s_obs + (n-1)/ n * (f1*(f1-1))/2
    }
  }
  return(sqrt((mean(s_chao)-S)^2 + var(s_chao)))
}
simu_func_proposed_RMSE <- function(n){
  thres = 0.9
  s_est <- 0
  for(i in c(1:500)){
    data <- rmultinom(1, n, pi_norm)
    sample_data <- data[data > 0]
    p_i <- sample_data / sum(sample_data)
    f1 <- sum(sample_data == 1)
    f2 <- sum(sample_data == 2)
    C <- 1-f1/n
    s_obs <- sum(sample_data > 0)
    if(f2 != 0){
      if(C < thres){
        s_est[i] <- sum(1/(1-(1-p_i)^n)) + (n-1)/ n * (f1^2)/(2*f2)
      }else{
        s_est[i] <- sum(1/(1-(1-p_i)^n))
      }
    }else{
      if(C < thres){
        s_est[i] <- sum(1/(1-(1-p_i)^n)) + (n-1)/ n * (f1*(f1-1))/2
      }else{
        s_est[i] <- sum(1/(1-(1-p_i)^n))
      }
    }
  }
  return(sqrt((mean(s_est)-S)^2 + var(s_est)))
}
vsimu_func_Chao <- Vectorize(simu_func_Chao)
vsimu_func_proposed <- Vectorize(simu_func_proposed)
vsimu_func_Chao_RMSE <- Vectorize(simu_func_Chao_RMSE)
vsimu_func_proposed_RMSE <- Vectorize(simu_func_proposed_RMSE)

plot(100:800, vsimu_func_Chao(100:800), col ="blue", type="l", ylim=c(100,250), 
     main = paste("Comparison of Species Richness Estimator", model, "CV:", round(cv,3)), xlab = "n", ylab = "Estimate")
lines(100:800, vsimu_func_proposed(100:800), type="l", col="red")
abline(h=S)

legend("bottomright", inset=.05, title="Estimator",
       c("Chao1","Proposed"), fill=c("blue", "red"), horiz=TRUE)

plot(100:800, vsimu_func_Chao_RMSE(100:800), col ="blue", type="l", ylim=c(0,100), 
     main = paste("Comparison of Species Richness Estimator", model, "CV:", round(cv,3)), xlab = "n", ylab = "RMSE")
lines(100:800, vsimu_func_proposed_RMSE(100:800), type="l", col="red")
abline(h=0)
legend("topright", inset=.05, title="Estimator",
       c("Chao1","Proposed"), fill=c("blue", "red"), horiz=TRUE)
print(data.frame(Model = model, CV = cv, True = S))
