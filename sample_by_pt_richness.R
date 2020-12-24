S <- 100
model <- 'unif'
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
n <- 50
x <- rmultinom(1, n, pi_norm)
chao <- numeric(0)
sobs <- numeric(0)
f1_s <- numeric(0)
f2_s <- numeric(0)
for(j in c(1:500)){
  n <- n + 1
  x <- x + rmultinom(1, 1, pi_norm)
  sobs[j] <- sum(x > 0)
  f1 <- sum(x == 1)
  f1_s[j] <- f1
  f2 <- sum(x == 2)
  f2_s[j] <- f2
  if(f2 == 0){
    chao[j] <- sum(x > 0) + (n-1)/(2*n) * f1 * (f1-1)
  }else{
    chao[j] <- sum(x > 0) + (n-1)/n * f1^2 / (2*f2)
  }
}
mean(f1_s[which(chao > S-1 & chao < S+1)]/f2_s[which(chao > S-1 & chao < S+1)])
plot(chao,type="l")
abline(h=S)
