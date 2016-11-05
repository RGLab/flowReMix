# Artificial data to work with ------------
n <- 100
nConditions <- 2
N <- sample(rv144$parentcount, n * nConditions)
ptid <- as.vector(replicate(n, rep(runif(1), nConditions)))
condition <- as.vector(replicate(n, 1:nConditions - 1))
randomsig <- 0.39
randomEffect <- as.vector(replicate(n, rep(rnorm(1, sd = randomsig), 2)))
M <- 10^3*5
simdata <- data.frame(ptid, condition, N)
X <- model.matrix(~ factor(condition), simdata)
beta <- c(-6, 0.7)
eta <- as.vector(X %*% beta) + randomEffect
expit <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log(x / (1 - x))
mu <- expit(eta)
p <- rbeta(length(mu), M*mu, M*(1 - mu))
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
  logdens <- lfactorial(N) - lfactorial(x) -
    lfactorial(N - x) - lbeta(beta, alpha) +
    lbeta(N - x + beta, x + alpha)
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

iters <- 1000
betaEst <- initBeta
Mest <- unique(imputedData$M[1])
randomsd <- numeric(iters)
MestVec <- numeric(iters)
for(i in 1:iters) {
  imputedData$eta[notNAind] <- X %*% betaEst

  randomSamp <- rnorm(400, sd = randomsig)
  randomEffectEstIter <- unlist(by(imputedData, ptid, computeRandomExpectation, randomSamp, simplify = TRUE))
  imputedData$randomEffectEst <- with(imputedData, randomEffectEst * .9 + 0.1 * randomEffectEstIter)

  imputedData$estEta <- with(imputedData, eta + randomEffectEst)
  imputedData$mu <- with(imputedData, expit(estEta))
  gradient <- with(imputedData, lbetabinomDerivMu(y, N, mu, M) * expitDeriv(estEta) * weights)
  gradient <- apply(cbind(gradient[notNAind], X), 1,
                    function(x) x[1] * x[2:(ncolX + 1)])
  gradient <- apply(gradient, 1, psych::winsor.mean) * ncol(gradient)

  betaEst <- betaEst - sign(gradient) * pmin(0.05, abs(gradient)) / sqrt(i)
  randomsd[i] <- with(imputedData, sd(randomEffectEst))
  #randomsig <- randomsig + 0.01/sqrt(i) * with(imputedData, ifelse(randomsig < psych::winsor.sd(randomEffectEst), 1, -1))
  randomsig <- randomsig + 0.01/sqrt(i) * with(imputedData, ifelse(randomsig < sqrt(mean(randomEffectEst^2)), 1, -1))
  print(sqrt(mean(imputedData$randomEffectEst^2)))

  imputedData$variance <- with(imputedData, mu * (1 - mu))
  imputedData$presiduals <- with(imputedData, (prop - mu) / sqrt(variance))
  Mest <- min(0.1, with(imputedData, weighted.mean(presiduals^2, weights, na.rm = TRUE)))
  Mest <- (1 - Mest)/Mest
  MestVec[i] <- Mest

  imputedData$M <- with(imputedData, Mest * 0.05 + M * .95)

  print(as.numeric(c(i, betaEst, imputedData$M[1], randomsig)))
}

with(imputedData, plot(trueRandomEffect, randomEffectEst, ylim = c(- 1, 1)))
abline(a = 0, b = 1)
plot(density(imputedData$randomEffectEst))
plot(density(imputedData$trueRandomEffect[notNAind]))


print(gradient)
print(psych::winsor.sd(imputedData$randomEffectEst))
