library(ggplot2)

S <- 200
r <- 1000
thres <- 1.5
model <- 'unif'
set.seed(1111)

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


for(model in c("unif", "zipf", "broken", "homo")){
  if(model == 'unif'){
    pi_i <- runif(S)
  }else if(model == 'zipf'){
    pi_i <- 1/((1:S)-0.1)
  }else if(model == 'broken'){
    pi_i <- rexp(S)
  }else if(model == 'homo'){
    pi_i <- rep(1/S, S)
  }
  pi_norm <- pi_i/sum(pi_i)
  cv <- (sd(pi_norm) / mean(pi_norm))
  
  richness <- ggplot(data.frame(x = c(100, 800)), aes(x)) +
    stat_function(fun = vsimu_func_Chao, aes(colour = "Chao1"), size=1.5) +
    stat_function(fun = vsimu_func_proposed, aes(colour = "Proposed"), size=1.5) +
    ylim(100, 280) +
    geom_hline(yintercept=S, color = "black") +
    ggtitle(paste("Model:",model, "; CV:", round(cv,3), "; S:", S)) +
    xlab("sample size") + 
    ylab("Estimation")
  print(richness)
  ggsave(filename = paste(model, "_richness.png", sep = ""), richness, dpi = 300, device='png')
  
  rmse <- ggplot(data.frame(x = c(100, 800)), aes(x)) +
    stat_function(fun = vsimu_func_Chao_RMSE, aes(colour = "Chao1"), size=1.5) +
    stat_function(fun = vsimu_func_proposed_RMSE, aes(colour = "Proposed"), size=1.5) +
    ylim(0, 100) +
    geom_hline(yintercept=0, color = "black") +
    ggtitle(paste("Model:", model, "; CV:", round(cv,3), "; S:", S)) +
    xlab("sample size") + 
    ylab("RMSE")
  print(rmse)
  ggsave(filename = paste(model, "_RMSE.png", sep = ""), rmse, dpi = 300, device='png')
  #print(data.frame(Model = model, CV = cv, True = S))
}

