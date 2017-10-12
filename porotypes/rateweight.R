iterations <- 1000
weights <- matrix(0, nrow = iterations, ncol = iterations)
weights[1, 1] <- 1
rate <- 0.55
for(i in 2:iterations) {
  weights[i, 1:(i - 1)] <- weights[i - 1, 1:(i-1)]  - weights[i - 1, 1:(i-1)] * 1/i^rate
  weights[i, i] <- 1/i^rate
}

sums <- t(apply(weights, 1, function(x) cumsum(x[iterations:1])))
discardAt <- 0.9
nsamps <- apply(sums, 1, function(x) sum(x != 0 & x < 0.9))
