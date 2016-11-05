# Artificial data to work with ------------
n <- 25
nConditions <- 2
N <- sample(rv144$parentcount, n * nConditions) * 100
ptid <- as.vector(replicate(n, rep(runif(1), nConditions)))
condition <- as.vector(replicate(n, 1:nConditions - 1))
trueRandomCov <- matrix(.95, nrow = nConditions, ncol = nConditions)
diag(trueRandomCov) <- 1
trueRandomCov <- trueRandomCov * 0.25
randomEffect <- as.vector(t(mvtnorm::rmvnorm(n, sigma = trueRandomCov)))
simdata <- data.frame(ptid, condition, N)
X <- model.matrix(~ factor(condition), simdata)
beta <- c(-6, 1)
eta <- as.vector(X %*% beta) + randomEffect
expit <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log(x / (1 - x))
mu <- expit(eta)
p <- mu
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
randomEffectMat <- with(imputedData, matrix(randomEffectEst, ncol = nWaves, byrow = TRUE))
randomCov <- cov(randomEffectMat, use = "pairwise.complete.obs")
imputedData$y <- imputedData$count
imputedData$normWeights <- with(imputedData, weights / sum(weights, na.rm = TRUE))

# Optimization loop
notNAind <- !is.na(imputedData$eta)
ncolX <- ncol(X)
MHsampler <- function(dataobj, randomCov, nSamp = 100, sampCoef = 0.4) {
  v <- dataobj$randomEffectEst
  prevRandom <- v
  y <- dataobj$y
  N <- dataobj$N
  eta <- dataobj$eta
  mu <- expit(eta + v)
  if(any(is.na(y))) {
    c("hello")
    notNA <- !is.na(prevRandom)
    randomCov <- randomCov[notNA, notNA, drop = FALSE]
    prevRandom <- prevRandom[notNA]
    y <- y[notNA]
    N <- N[notNA]
    eta <- eta[notNA]
    mu <- mu[notNA]
    v <- v[notNA]
  } else {
    notNA <- rep(TRUE, length(v))
  }
  prevlogdens <- sum(dbinom(y, N, mu, log = TRUE), na.rm = TRUE) + mvtnorm::dmvnorm(v, sigma = randomCov, log = TRUE)
  randomEffectSamp <- matrix(ncol = nrow(randomCov), nrow = nSamp)
  vCandidates <- mvtnorm::rmvnorm(nSamp, mean = prevRandom, sigma = randomCov * sampCoef)
  for(i in 1:nSamp) {
    v <- vCandidates[i, ] + prevRandom
    mu <- expit(eta + v)

    logdens <- sum(dbinom(y, N, mu, log = TRUE), na.rm = TRUE) + mvtnorm::dmvnorm(v, sigma = randomCov, log = TRUE)
    MHratio <- exp(logdens - prevlogdens)
    #print(rbind(unique(dataobj$trueRandomEffect), prevRandom, v))
    #print(MHratio)
    U <- runif(1)
    test <- U < MHratio
    try(if(test) {
      prevRandom <- v
      prevlogdens <- logdens
    }, silent = TRUE)
    randomEffectSamp[i, ] <- prevRandom
  }
  #print(colMeans(randomEffectSamp))
  #print(with(dataobj, trueRandomEffect))
  #return(rep(colMeans(randomEffectSamp), length(y)))

  result <- numeric(length(notNA))
  result[!notNA] <- NA
  result[notNA] <- colMeans(randomEffectSamp)
  print(cbind(y/N, expit(eta + result)))
  return(result)
}

nIters <- 200
nSamp <- 10
notNAind <- !is.na(imputedData$y)
for(i in 1:nIters) {
  #system.time(MHresult <- by(imputedData, ptid, MHsampler, randomsig = 0.4, nSamp = nSamp, simplify = TRUE))
  #MHresult <- matrix(unlist(MHresult), byrow = TRUE, ncol = nSamp)
  #St <- rowMeans(MHresult)
  system.time(St <- unlist(by(imputedData, imputedData$ptid, MHsampler,
                              randomCov = randomCov,
                              nSamp = nSamp,
                              sampCoef = 0.2,
                              simplify = TRUE)))
  newRandomEst <- with(imputedData, randomEffectEst +  (St - randomEffectEst) / sqrt(i))
  imputedData$randomEffectEst <- newRandomEst

  system.time(regfit <- glm(prop ~ X - 1 + offset(randomEffectEst), data = subset(imputedData, notNAind),
                family = "quasibinomial"))
  imputedData$eta[notNAind] <- predict(regfit)
  beta <- coef(regfit)
  newCov <- cov(matrix(St, ncol = nWaves, byrow = TRUE), use = "pairwise.complete.obs")
  randomCov <- randomCov + 0.5 * (newCov - randomCov) / sqrt(i)

  print(i)
  print(beta)
  print((randomCov))
  with(imputedData, print(mean((randomEffectEst - trueRandomEffect)^2, na.rm = TRUE)))
  #apply(MHresult, 1, function(x) mean((x[2:nSamp] - x[1:(nSamp - 1)]) != 0 ))
}

par(mfrow = c(2, 1), mar = rep(3, 4))
with(imputedData, plot(density(randomEffectEst[notNAind])))
range <- seq(from = -3 * trueRandomCov[1, 1], to = 3 * trueRandomCov[1, 1], length.out = 10^3)
lines(range, dnorm(range, sd = sqrt(trueRandomCov[1, 1])), col = "red")
with(imputedData, plot(randomEffectEst, trueRandomEffect))
abline(a = 0, b = 1)
with(imputedData, summary(lm(randomEffectEst ~ trueRandomEffect)))
