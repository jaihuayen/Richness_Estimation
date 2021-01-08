S <- 200
r <- 200
model <- 'zipf'
set.seed(1111)

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

chao_fun <- function(data){
  sample_data <- data[data > 0]
  s_obs <- sum(sample_data > 0)
  f1 <- sum(sample_data == 1)
  f2 <- sum(sample_data == 2)
  if(f2 > 0){
    s_chao <- s_obs + (n-1)/ n * (f1^2)/(2*f2)
  }else{
    s_chao <- s_obs + (n-1)/(2*n) * f1 * (f1-1)
  }
  return(s_chao)
}

jk_chao_fun <- function(data){
  sample_data <- data[data > 0]
  jk_chao <- 0
  for(i in 1:length(sample_data)){
    sample_data_jk <- sample_data[-i]
    jk_chao[i] <- chao_fun(sample_data_jk)
  }
  return(mean(jk_chao))
}

simu_func_Chao <- function(n){
  chao <- 0
  for(i in c(1:r)){
    data <- rmultinom(1, n, pi_norm)
    chao[i] <- chao_fun(data)
  }
  return(mean(chao))
}
vsimu_func_Chao <- Vectorize(simu_func_Chao)

simu_func_jk_Chao <- function(n){
  jk_chao <- 0
  for(i in c(1:r)){
    data <- rmultinom(1, n, pi_norm)
    jk_chao[i] <- jk_chao_fun(data)
  }
  return(mean(jk_chao))
}
vsimu_func_jk_Chao <- Vectorize(simu_func_jk_Chao)

sample_size_vec <- 50:800

plot(sample_size_vec,vsimu_func_Chao(c(sample_size_vec)), type = 'l', 
     col='black', ylim = c(S*0.5,S*1.5), lwd = 4,
     main = paste("Model:",model, "; CV:", round(cv,3), "; S:", S))
lines(sample_size_vec,vsimu_func_jk_Chao(c(sample_size_vec)), type = 'l', 
      col='red', lwd = 4)
abline(h=S,col="dark green",lty=2)
