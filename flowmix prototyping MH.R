# Artificial data to work with ------------
n <- 100
nConditions <- 2
N <- sample(rv144$parentcount, n * nConditions) * 10
ptid <- as.vector(replicate(n, rep(runif(1), nConditions)))
condition <- as.vector(replicate(n, 1:nConditions - 1))
randomsig <- 0.1
trueRandomsig <- randomsig
randomEffect <- as.vector(replicate(n, rep(rnorm(1, sd = randomsig), nConditions)))
M <- 10^4
trueM <- M
simdata <- data.frame(ptid, condition, N)
X <- model.matrix(~ factor(condition), simdata)
beta <- c(-5, 0.5)
eta <- as.vector(X %*% beta) + randomEffect
expit <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log(x / (1 - x))
mu <- expit(eta)
p <- rbeta(length(mu), M*mu, M*(1 - mu))
boxplot(p ~ factor(expit(as.vector(X %*% beta))))
y <- rbinom(length(mu), N, p)
logpropmat <- matrix((y/N), ncol = nConditions, byrow = TRUE)
plot(logpropmat[, 1], logpropmat[, 2],
     ylim = c(min(logpropmat), max(logpropmat)),
     xlim = c(min(logpropmat), max(logpropmat)))
abline(a = 0, b = 1)
simdata$count <- y
waves <- as.vector(replicate(n, 1:nConditions))
simdata$waves <- waves
simdata$trueRandomEffect <- randomEffect
simdata <- simdata[-1, ]
simdata$prop <- simdata$count / simdata$N
simdata$weights <- rep(1, nrow(simdata))

#Pre-processing
imputedData <- with(simdata, expand.grid(waves = unique(waves), ptid = unique(ptid)))
imputedData <- merge(imputedData, simdata, all = TRUE)
dataOrder <- with(imputedData, order(ptid, waves))
imputedData <- imputedData[dataOrder, ]
X <- model.matrix(~ factor(condition), imputedData)
nWaves <- length(unique(waves))
imputedData$prop <- imputedData$count / imputedData$N
notMissingInd <- !is.na(imputedData[, 3])

# Initial Estimate for Beta
initFit <- glm(cbind(count, N - count) ~ X - 1,
              family = "binomial",
              data = subset(imputedData, notMissingInd),
              weights = weights)
initBeta <- coef(initFit)

# Initial Random Effect Estimate
imputedData$empEta <- with(imputedData, logit(prop))
imputedData$eta[notMissingInd] <- predict(initFit)
imputedData$randomEffectEst <- with(imputedData, eta - empEta)
etaDiffMat <- with(imputedData,
                   matrix(eta - empEta, byrow = TRUE, ncol = nWaves))
randomEffectEst <- rowMeans(etaDiffMat, na.rm = TRUE)
randomsig <- sd(randomEffectEst)
initRandomsig <- randomsig
randomEffectEst <- as.vector(sapply(randomEffectEst, function(x) rep(x, nWaves)))
imputedData$randomEffectEst <- rowMeans(etaDiffMat, na.rm = TRUE)
#imputedData$mu <- with(imputedData, expit(randomEffectEst + eta))
imputedData$mu <- with(imputedData, expit(eta))
imputedData$variance <- with(imputedData, mu * (1 - mu))
imputedData$residuals <- with(imputedData, prop - mu)
imputedData$presiduals <- with(imputedData, residuals / sqrt(variance))
rho <- with(imputedData, weighted.mean(presiduals^2, weights, na.rm = TRUE))
imputedData$M[notMissingInd] <- (1 - rho) / rho
imputedData$y <- imputedData$count
imputedData$normWeights <- with(imputedData, weights / sum(weights, na.rm = TRUE))

# Some functions for optimization --------------------
dbetabinom <- function(x, N, mu, M, log = TRUE) {
  alpha <- M * mu
  beta <- M * (1 - mu)
  logdens <- lfactorial(N) - lfactorial(x) - lfactorial(N - x) -
    lbeta(M * (1 - mu), M * mu) +
    lbeta(N - x + M * (1 - mu), x + M * mu)
  if(!log) return(exp(logdens))
  return(logdens)
}
lbetabinomDerivMu <- function(x, N, mu, M) {
  MtimesMu <- M * mu
  deriv <- -digamma(M * mu) * M
  deriv <- deriv + M * digamma(M - MtimesMu)
  deriv <- deriv + digamma(x + MtimesMu) * M
  deriv <- deriv - M * digamma(N - x + M - MtimesMu)
  return(deriv)
}
lbetabinomDerivDispersion <- function(x, N, mu, M) {
  oneMinusMu <- (1 - mu)
  MtimesMu <- M * mu
  MtimesOneMinusMu <- M * oneMinusMu
  deriv <- digamma(M)
  deriv <- deriv - mu * digamma(MtimesMu)
  deriv <- deriv - digamma(MtimesOneMinusMu) * oneMinusMu
  deriv <- deriv + digamma(x + MtimesMu) * mu
  deriv <- deriv + digamma(N - x + MtimesOneMinusMu) * oneMinusMu
  deriv <- deriv - digamma(N + M)
  return(deriv)
}
expitDeriv <- function(eta) {
  deriv <- 1/(1 + exp(- eta))^2
}

# Optimization loop
notNAind <- !is.na(imputedData$eta)
ncolX <- ncol(X)
computeRandomExpectation <- function(dataobj, randomSamp) {
  eta <- dataobj$eta
  M <- dataobj$M
  N <- dataobj$N
  x <- dataobj$y
  sampSize <- length(randomSamp)
  likelihood <- sapply(randomSamp, function(v)
    sum(dbetabinom(x,  N, expit(eta + v),  M, log = TRUE), na.rm = TRUE))
  likelihood <- likelihood - max(likelihood)
  denominator <- exp(likelihood[1:round(sampSize / 2)])
  numeratorIndex <- (round(sampSize/2) + 1):sampSize
  numerator <- (exp(likelihood[numeratorIndex]) * randomSamp[numeratorIndex])
  expNum <- mean(numerator)
  expDenom <- mean(denominator)
  varDenom <- var(denominator)
  ratioExpectation <- expNum / expDenom - varDenom * expNum / expDenom^3
  return(rep(ratioExpectation, nrow(dataobj)))
}

MHsampler <- function(dataobj, randomsig, nSamp = 100, sampsig = NULL) {
  if(is.null(sampsig)) sampsig <- randomsig
  v <- dataobj$randomEffectEst[1]
  prevRandom <- v
  y <- dataobj$y
  N <- dataobj$N
  eta <- dataobj$eta
  M <- dataobj$M
  mu <- expit(eta + v)
  prevlogdens <- sum(dbetabinom(y, N, mu, M), na.rm = TRUE) + dnorm(v, sd = randomsig, log = TRUE)
  randomEffectSamp <- numeric(nSamp)
  for(i in 1:nSamp) {
    v <- rnorm(1, mean = prevRandom, sd = sampsig)
    mu <- expit(eta + v)

    logdens <- sum(dbetabinom(y, N, mu, M), na.rm = TRUE) + dnorm(v, sd = randomsig, log = TRUE)
    MHratio <- exp(logdens - prevlogdens)
    #print(c(unique(dataobj$trueRandomEffect), prevRandom, v, MHratio))
    #print(cbind(y/N, mu))
    U <- runif(1)
    test <- U < MHratio
    try(if(test) {
      prevRandom <- v
      prevlogdens <- logdens
    }, silent = TRUE)
    randomEffectSamp[i] <- prevRandom
    #print(round(c(unique(dataobj$trueRandomEffect), prevRandom, v, MHratio, U, prevlogdens), 4))
  }
  #print(cbind(y/N, expit(eta + prevRandom)))
  #print(mean(randomEffectSamp))
  #return(rep(mean(randomEffectSamp), length(y)))
  return(rep(prevRandom, length(y)))
}

nIters <- 50
nSamp <- 40
notNAind <- !is.na(imputedData$y)
rhoest <- 10^-4
for(i in 1:nIters) {
  #system.time(MHresult <- by(imputedData, ptid, MHsampler, randomsig = 0.4, nSamp = nSamp, simplify = TRUE))
  #MHresult <- matrix(unlist(MHresult), byrow = TRUE, ncol = nSamp)
  #St <- rowMeans(MHresult)
  system.time(St <- unlist(by(imputedData, imputedData$ptid, MHsampler,
                              randomsig = randomsig * 1,
                              nSamp = nSamp,
                              sampsig = trueRandomsig,
                              simplify = TRUE)))
  #St <- randomEffect
  newRandomEst <- with(imputedData, randomEffectEst +  (St - randomEffectEst) / i)
  #newRandomEst <- randomEffect
  imputedData$randomEffectEst <- newRandomEst

  system.time(regfit <- glm(prop ~ X - 1 + offset(randomEffectEst), data = subset(imputedData, notNAind),
                family = "quasibinomial"))
  imputedData$eta[notNAind] <- predict(regfit)
  beta <- coef(regfit)
  randomsig <- randomsig + 0.5 * (sd(St) - randomsig) / i
  rho <- summary(regfit)$dispersion
  M <- (1 - rho)/rho
  imputedData$M <- M

  muest <- predict(regfit, type = "response")
  varest <- muest*(1-muest)
  presid <- with(subset(imputedData, notNAind), (prop - muest)/sqrt(varest))
  rhoest <- mean(presid^2)
  binomvar <- varest / N[notNAind]
  betavar <- varest
  rhoest <- mean(presid^2)
  Mest <- (1 - rhoest)/rhoest
  imputedData$M <- Mest

  print(i)
  print(beta)
  print(c(randomsig, rho, 1/(1 + trueM)))
  with(imputedData, print(mean((randomEffectEst - trueRandomEffect)^2, na.rm = TRUE)))
  #apply(MHresult, 1, function(x) mean((x[2:nSamp] - x[1:(nSamp - 1)]) != 0 ))
}

par(mfrow = c(2, 1), mar = rep(3, 4))
with(imputedData, plot(density(randomEffectEst)))
range <- seq(from = -3 * trueRandomsig, to = 3*trueRandomsig, length.out = 10^3)
lines(range, dnorm(range, sd = trueRandomsig), col = "red")
with(imputedData, plot(randomEffectEst, trueRandomEffect))
abline(a = 0, b = 1)
with(imputedData, summary(lm(randomEffectEst ~ trueRandomEffect)))
