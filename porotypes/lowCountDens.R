betaBinomDens <- function(y, N, M, prob) {
  a <- M * prob ;
  b <- M * (1 - prob) ;
  logdens <- lchoose(N, y) + lbeta(y + a, N - y + b) - lbeta(a, b)
  return(logdens)
}


load(file = "data analysis/results/boolean robust5 wrandom.Robj")
subset <- 23
prob0 <- 1 / (1 + exp(-fit$coefficients[[subset]][1]))
prob1 <- 1 / (1 + exp(-sum(fit$coefficients[[subset]][1:2])))
N <- 10^4
y0 <- c(0, 0, 0, 0)
y1 <- c(0, 0, 1, 3)
M <- 10^9
prob <- 1
range <- 0:max(1, N * prob1)
dens0 <- exp(betaBinomDens(range, N, M, prob0))
dens1 <- exp(betaBinomDens(range, N, M, prob1))
plot(range, dens0, type = "l", ylim = c(min(dens0, dens1), max(dens0, dens1)))
lines(range, dens1, type = "l", col = "red")
cbind(dens0, dens1)
posprob <- 0.6
dens0[1] * dens0[1] * (1 - posprob)
dens0[1] * dens1[1] * posprob

dbinom(0, 10^4, prob0)
dbinom(0, 10^4, prob1)
