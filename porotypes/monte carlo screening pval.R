reps <- 1000
ratios <- numeric(reps)
nullprob <- 0.5
M <- 100
n <- 260
N <- 600
N <- rep(N, n)

mcScreenTest <- function(y, N, reps = 1000) {
  n <- length(y)
  if(length(N) == 1) {
    N <- rep(N, n)
  } else {
    if(length(N) != n) stop("wrong length for N!")
  }

  nullprop <- rep(mean(y / N), n)
  M <- max(dispersionMLE(y, N, nullprop), 1000)
  props <- y / N
  nullprop <- rep(mean(props), n)
  ratio <- sum(dbinom(y, N, props, TRUE)) - sum(vecBetaBinomDens(y, N, nullprop, M))

  ratios <- numeric(reps)
  for(i in 1:reps) {
    probs <- rbeta(n, nullprob * M, (1 - nullprob) * M)
    y <- rbinom(n, N, probs)
    props <- y/N
    nullprop <- rep(mean(y / N), n)
    estM <- dispersionMLE(y, N, nullprop)
    ratios[i] <- sum(dbinom(y, N, props, TRUE)) - sum(vecBetaBinomDens(y, N, nullprop, estM))
  }

  pval <- 1- mean(ratio > ratios)
  return(list(pval = pval, M = M))
}
