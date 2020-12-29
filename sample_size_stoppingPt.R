S <- 200
t <- 100
model <- 'homo'

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

criterion_function <- function(n, f1, f2){
  if(f2 == 0){
    d <- (n-1)/(2*n)*f1*(f1-1)*(1-2/((n-1)*(f1-1)+2))
    return(d)
  }else{
    d <- (n-1)/(2*n*f2)*(f1^2)*(1-2*f2/((n-1)*f1+f2))
    return(d)
  }
}


nn <- chao <- numeric(0)
for(j in c(1:t)){
  n <- 50
  x <- rmultinom(1, n, pi_norm)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  while(criterion_function(n, f1, f2) >= 1){
    n <- n + 1
    x <- x + rmultinom(1, 1, pi_norm)
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
  }
  if(f2 == 0){
    chao[j] <- sum(x > 0) + (n-1)/(2*n) * f1 * (f1-1)
  }else{
    chao[j] <- sum(x > 0) + (n-1)/n * f1^2 / (2*f2)
  }
  
  nn[j] <- n
}
app1 <- data.frame(Method = "Rarefaction",n = mean(nn), est = mean(chao),
                   se_n = sd(nn), rmse = sqrt((S-mean(chao))^2+var(chao)))



nn <- chao <- numeric(0)
for(j in c(1:t)){
  n <- 50
  x <- rmultinom(1, n, pi_norm)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  while(f1 > 0){
    n <- n + 1
    x <- x + rmultinom(1, 1, pi_norm)
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
  }
  if(f2 == 0){
    chao[j] <- sum(x > 0) + (n-1)/(2*n) * f1 * (f1-1)
  }else{
    chao[j] <- sum(x > 0) + (n-1)/n * f1^2 / (2*f2)
  }
  nn[j] <- n
}
app2 <- data.frame(Method = "f0=0",n = mean(nn), est = mean(chao),
                   se_n = sd(nn), rmse = sqrt((S-mean(chao))^2+var(chao)))
nn <- chao <- numeric(0)
for(j in c(1:t)){
  n <- 50
  x <- rmultinom(1, n, pi_norm)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  f0_fun <- function(n,f1,f2){
    if(f2 == 0){
      return((n-1)/(2*n) * f1 * (f1-1))
    }else{
      return((n-1)/n * f1^2 / (2*f2))
    }
  }
  while(f0_fun(n,f1,f2) >= 1){
    n <- n + 1
    x <- x + rmultinom(1, 1, pi_norm)
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
  }
  if(f2 == 0){
    chao[j] <- sum(x > 0) + (n-1)/(2*n) * f1 * (f1-1)
  }else{
    chao[j] <- sum(x > 0) + (n-1)/n * f1^2 / (2*f2)
  }
  nn[j] <- n
}
app3 <- data.frame(Method = "f0<1",n = mean(nn), est = mean(chao),
                   se_n = sd(nn), rmse = sqrt((S-mean(chao))^2+var(chao)))

nn <- chao <- numeric(0)
for(j in c(1:t)){
  n <- 50
  x <- rmultinom(1, n, pi_norm)
  sobs <- sum(x > 0)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  if(f2 == 0){
    while((sobs - (n-1)/n * (f1 * (f1-1) / 2)) <= 0){
      n <- n + 1
      x <- x + rmultinom(1, 1, pi_norm)
      sobs <- sum(x > 0)
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
    }
  }else{
    while((sobs - (n-1)/n *(f1^2 / 2*f2)) <= 0){
      n <- n + 1
      x <- x + rmultinom(1, 1, pi_norm)
      sobs <- sum(x > 0)
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
    }
  }

  if(f2 == 0){
    chao[j] <- sum(x > 0) + (n-1)/(2*n) * f1 * (f1-1)
  }else{
    chao[j] <- sum(x > 0) + (n-1)/n * f1^2 / (2*f2)
  }
  nn[j] <- n
}
app4 <- data.frame(Method = "sobs > f1^2 / 2*f2",n = mean(nn), est = mean(chao),
                   se_n = sd(nn), rmse = sqrt((S-mean(chao))^2+var(chao)))

rbind(app1,app2,app3,app4)
