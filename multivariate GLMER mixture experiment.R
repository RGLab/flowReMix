library(flowReMix)
cummean <- function(x) cumsum(x) / 1:length(x)
data(rv144)
#set.seed(502)
set.seed(504)
par(mfrow = c(1, 1), mar = rep(4, 4))
data <- rv144
leaves <- unique(data$population)
data <- subset(data, population %in% leaves[c(1:3, 5:7)])
data <- subset(data, stim != "sebctrl")
data$treatment <- as.numeric(data$stim == "env")
data$ptid <- as.numeric(data$ptid)
data$ptid[data$vaccine == "VACCINE"] <- data$ptid[data$vaccine == "VACCINE"] * 10^4
data$prop <- data$count / data$parentcount
data$population <- as.factor(data$population)
data <- data[order(data$population, data$ptid, data$stim, decreasing = FALSE), ]

mixtureFitList <- by(data, data$population, function(X)
                     mixedfit <- glmmMixture(count ~ (age + gender + treatment),
                                             N = parentcount,
                                             sub.population = population,
                                             id = ptid,
                                             treatment = treatment,
                                             data = X,
                                             tol = 0.01,
                                             maxiter = 1,
                                             nAGQ = 1))
# Getting list of Coefficients
coefficientList <- lapply(mixtureFitList, function(x) x$beta)

# Estimating covariance structure from marginal fits (step 0)
randomEffects <- lapply(mixtureFitList, function(x) as.vector(x$randomEffectEst))
weights <- lapply(mixtureFitList, function(x) (x$subject.probs))
levelProbs <- colMeans(do.call("rbind", weights))
weights <- lapply(weights, function(x) as.vector(x))
weights <- do.call("cbind", weights)
weights <- rowMeans(weights)
randomEffects <- do.call("cbind", randomEffects)
covariance <- cov.wt(randomEffects, weights, center = FALSE)$cov
round(covariance, 3)
round(cov2cor(covariance), 3)

# Computing new posterior porbabilities
maxIter <- 30
vaccines <- sapply(databyid, function(x) x$vaccine[1] == "VACCINE")
muMat <- lapply(mixtureFitList, function(x) x$mu[, 4:5])
muMat <- do.call("rbind", muMat)
data$nullMu <- muMat$nullMu
data$altMu <- muMat$altMu
data$subpopInd <- as.numeric(data$population)
databyid <- by(data, data$ptid, function(x) x)
nSubjects <- length(unique(data$ptid))
sampCoef <- 0.00001
sampcov <- sampCoef * covariance
sqrtcov <- expm::sqrtm(sampcov)
invcov <- solve(covariance)
invSampcov <- solve(sampcov)
covDet <- as.numeric(determinant(covariance, logarithm = TRUE)$modulus)
sampCovDet <- as.numeric(determinant(sampcov, logarithm = TRUE)$modulus)
allPosteriors <- matrix(nrow = nSubjects, ncol = maxIter)
estimatedRandomEffects <- randomEffects
randomSampList <- lapply(1:2, function(x) x)
posteriors <- numeric(nSubjects)
clusterAssignments <- numeric(nSubjects)
rate <- 1
lastMean <- randomEffects
iterCoefMat <- matrix(ncol = length(mixtureFitList), nrow = maxIter + 1)
accept <- 0
for(iter in 1:maxIter) {
  nsamp <- 100 + iter
  logLikelihoods <- matrix(nrow = 2, ncol = nsamp)
  zSamp <- matrix(rnorm(nsamp * ncol(randomEffects)), nrow = ncol(randomEffects))
  randomEffectSamp <- sqrtcov %*% zSamp
  accept = flowReMix:::zero(accept)
  for(i in 1:nSubjects) {
    # Computing Posterior Probabilities (Step 1a)
    subjectData <- databyid[[i]]
    popInd <- subjectData$subpopInd
    N <- subjectData$parentcount
    y <- subjectData$count
    prop <- y/N
    for(k in 1:2) {
      randEst <- estimatedRandomEffects[2*i - 2 + k, ]
      randSamp <- randomEffectSamp+randEst #faster Sampapply(randomEffectSamp, 2, function(x) x + randEst)
      if(iter == 1) {
        if(k == 1) {
          eta <- logit(subjectData$nullMu) - randEst[popInd]
          subjectData$nullEta <- eta
          nullEta <- eta
        } else {
          eta <- logit(subjectData$altMu) - randEst[popInd]
          subjectData$altEta <- eta
          altEta <- eta
        }
      } else {
        if(k == 1) {
          eta <- subjectData$nullEta
          nullEta <- eta
        } else {
          eta <- subjectData$altEta
          altEta <- eta
        }
      }
      mu <- expit(eta+randSamp[popInd,]) #faster than expit(apply(randSamp, 2, function(x) eta + x[popInd])) and equivalent
      binomLlik <- colSums(dbinom(y,N,mu,log=TRUE)) # Faster than this: apply(mu, 2, function(x) sum(dbinom(y, N, x, log = TRUE)))
      normWeights <- #apply(randSamp, 2, function(x) - 0.5 * t(x) %*% invcov  %*% x - 0.5 * covDet)
              -0.5*diag(crossprod(crossprod(invcov,randSamp),randSamp))-0.5*covDet #faster than above
      importanceWeights <- #apply(randSamp, 2, function(x) - 0.5 * t(x - randEst) %*% invSampcov  %*% (x - randEst) - 0.5 * sampCovDet)
              -0.5*diag(crossprod(crossprod(invSampcov,randSamp-randEst),randSamp-randEst))-0.5*sampCovDet #faster than above
      logLikelihoods[k, ] <- binomLlik + normWeights - importanceWeights + log(levelProbs[k])
      randomSampList[[k]] <- randSamp
    }

    posterior <- rowSums(exp(logLikelihoods - max(logLikelihoods)))
    posterior <-  1/(1 + posterior[1] / posterior[2])
    if(is.nan(posterior)) posterior <- 0
    #print(c(i, subjectData$vaccine[1] == "VACCINE", posterior))
    if(iter == 1) {
      posteriors[i] <- posterior
    } else if(iter > 1) {
      posteriors[i] <- posteriors[i] + (posterior - posteriors[i]) / (iter)^rate
    }

    # Sampling cluster assignment (Step 1b)
    cluster <- 1 + rbinom(1, 1, posteriors[i])
    clusterAssignments[i] <- cluster
    subjectData$tempTreatment <- subjectData$treatment * (cluster == 2)

    # Performing MH step
    # This loop *may* be faster in C but it's hard to tell.. there's not that much sampling going on and gains may be lost due to the overhead of calling out.
    # This is the loop below implemented directly in C++. No real optimization was done yet, just a direct translation.
    unifs = runif(nsamp)
    # browser()
    # if(iter==2){stop();}
    flowReMix:::MH(lastMean, estimatedRandomEffects, y, N, randomEffectSamp, i, popInd, invcov, accept, iter,rate, unifs, nullEta, altEta);
    # for(k in 1:2) {
    #   currentSamp <- lastMean[2*i - 2 + k, ]
    #   mu <- expit(eta + currentSamp[popInd])
    #   currentloglik <- sum(dbinom(y, N, mu, log = TRUE)) - 0.5 * t(currentSamp) %*% invcov %*% (currentSamp)
    #   for(j in 1:ncol(randomEffectSamp)) {
    #     newSamp <- randomEffectSamp[, j] + currentSamp
    #     mu <- expit(eta + newSamp[popInd])
    #     newlogLik <- sum(dbinom(y, N, mu, log = TRUE)) - 0.5 * t(newSamp) %*% invcov %*% (newSamp)
    #     if(unifs[j] < exp(newlogLik - currentloglik)) {
    #       currentSamp <- newSamp
    #       currentloglik <- newlogLik
    #       accept <- accept + 1
    #     }
    #   # cat("Accepted ",accept," on subject ",i,"\n");
    #   lastMean[2*i - 2 + k, ] <- currentSamp
    #   currentEst <- estimatedRandomEffects[2*i - 2 + k, ]
    #   estimatedRandomEffects[2*i - 2 + k, ] <- currentEst + (currentSamp - currentEst) / (iter + 1.0)^rate
    #   }
    # }

    subjectData$randomOffset <- lastMean[2*i - 2 + cluster, ]
    #subjectData$randomOffset <- estimatedRandomEffects[2*i - 2 + cluster, ]
    databyid[[i]] <- subjectData
  }
  # Refitting Model with current means
  dataForGlm <- data.frame(data.table::rbindlist(databyid)) # much faster than:  do.call("rbind", databyid)
  dataByPopulation <- by(dataForGlm, dataForGlm$population, function(x) x)
  glmFits <- lapply(dataByPopulation, function(popdata)
                glm(cbind(count, parentcount - count) ~  age + gender + tempTreatment + offset(randomOffset),
                    family = "binomial", data = popdata))
  coefficientList <- mapply(function(coef, fit) coef + (coef(fit) - coef)/(iter)^rate,
                            coefficientList, glmFits, SIMPLIFY = FALSE)
  # Updating Prediction
  for(j in 1:length(dataByPopulation)) {
    dataByPopulation[[j]]$tempTreatment <- 0
    dataByPopulation[[j]]$randomOffset <- 0
    newNullEta <- predict(glmFits[[j]], newdata = dataByPopulation[[j]])
    dataByPopulation[[j]]$tempTreatment <- dataByPopulation[[j]]$treatment
    newAltEta <- predict(glmFits[[j]], newdata = dataByPopulation[[j]])
    nullEta <- dataByPopulation[[j]]$nullEta
    altEta <- dataByPopulation[[j]]$altEta
    dataByPopulation[[j]]$nullEta <- nullEta + (newNullEta - nullEta)/(iter)^rate
    dataByPopulation[[j]]$altEta <- altEta + (newAltEta - altEta)/(iter)^rate
  }

  databyid <- do.call("rbind", dataByPopulation)
  databyid <- with(databyid, databyid[order(population, ptid, stim, decreasing = FALSE), ])
  databyid <- by(databyid, databyid$ptid, function(x) x)

  acceptRate <- accept / (2*nsamp * nSubjects)
  # cat("acceptance rate: ",acceptRate,"\n")
  # Updating covariance
  if(iter == 1) {
    sampCoef <- 0.1
    sampcov <- covariance * sampCoef
    sqrtcov <- expm::sqrtm(sampcov)
    invSampcov <- solve(sampcov)
    sampCovDet <- as.numeric(determinant(sampcov, logarithm = TRUE)$modulus)
  } else {
    if(acceptRate < .234) {
      sampCoef <- sampCoef * .96
    } else if (acceptRate > 0.234) {
      sampCoef <- sampCoef * 1.03
    }
    sampcov <- covariance * sampCoef
    sqrtcov <- expm::sqrtm(sampcov)
    invSampcov <- solve(sampcov)
    sampCovDet <- as.numeric(determinant(sampcov, logarithm = TRUE)$modulus)
  }

  weights <- as.vector(sapply(posteriors, function(x) c(1 - x, x)))
  covariance <- cov.wt(estimatedRandomEffects, weights, center = FALSE)$cov
  # print(covariance)
  invcov <- solve(covariance)
  covDet <- as.numeric(determinant(covariance, logarithm = TRUE)$modulus)

  levelProbs[2] <- mean(posteriors)
  levelProbs[1] <- 1 - levelProbs[2]
  print(c(iter, sampCoef, accept / (2*nsamp * nSubjects)))
  print(levelProbs)

  # Some Diagnostics/Outputs
  iterCoef <- sapply(glmFits, function(x) coef(x)[5])
  iterCoefMat[iter + 1, ] <- iterCoef
  currentCoef <- sapply(coefficientList, function(x) x[[5]])
  initCoef <- sapply(mixtureFitList, function(x) x$beta[[5]])
  iterCoefMat[1, ] <- initCoef
  print(round(cbind(iterCoef, currentCoef, initCoef), 3))
  # print(round(cov2cor(covariance), 3))
  require(pROC)
  rocfit <- roc(vaccines ~ posteriors)
  print(plot(rocfit, main = round(rocfit$auc, 3)))
}










require(pROC)
vaccines <- sapply(databyid, function(x) x$vaccine[1] == "VACCINE")
rocfit <- roc(vaccines ~ posteriors)
plot(rocfit)


require(ggplot2)
populations <- unique(data$population)
plotList <- lapply(1:length(mixtureFitList), function(x) x)
for(i in 1:length(mixtureFitList)) {
  tempdat <- subset(data, population == populations[i])
  prop <- log(tempdat$count / tempdat$parentcount)
  stimprop <- prop[tempdat$stim == "env"]
  ctrlprop <- prop[tempdat$stim != "env"]
  vaccine <- (tempdat$vaccine == "VACCINE")[tempdat$stim == "env"]
  plotList[[i]] <- data.frame(vaccine = vaccine, stimprop = stimprop,
                              ctrlprop = ctrlprop, population = populations[i],
                              posterior = posteriors)
}

forPlot <- do.call("rbind", plotList)
require(ggplot2)
ggplot(forPlot) +
  geom_point(aes(x = ctrlprop, y = stimprop, col = posterior, shape = !vaccine), fill = "white") +
  theme_bw() + geom_abline(intercept = 0, slope = 1) +
  scale_colour_gradientn(colours=rainbow(4)) +
  facet_wrap(~ population, scales = "free", ncol = 3)

par(mfrow = c(1, 2), mar = rep(3, 4))
sprobs <- data.frame(posteriors, vaccine)
sprobs <- sprobs[order(posteriors, decreasing = TRUE), ]
sprobs$nominalFDR <- cummean(1 - sprobs$posteriors)
sprobs$empFDR <- cummean(1 - sprobs$vaccine)
sprobs$power <- cumsum(sprobs$vaccine) / sum(sprobs$vaccine)
uniqueNominal <- unique(sprobs$nominalFDR)
empFDR <- sapply(uniqueNominal, function(x) sprobs$empFDR[max(which(sprobs$nominalFDR == x))])
power <- sapply(uniqueNominal, function(x) sprobs$power[max(which(sprobs$nominalFDR == x))])
lim <- max(c(empFDR, uniqueNominal))
plot(uniqueNominal, empFDR, type = "l", xlim = c(0, lim), ylim = c(0, 1), col = "red",
     xlab = "nominal FDR", ylab = "Empirical FDR / Power")
abline(a = 0, b = 1)
lines(uniqueNominal, power, col = "blue", lty = 2)
legend("topright", col = c("red", "blue"), lty = 1:2, legend = c("FDR", "Power"))
abline(v = c(.01, .05, .1), h = c(.75, .9, .95), col = "grey")
print(plot(rocfit, main = round(rocfit$auc, 3)))

require(xtable)
covTable <- xtable(covariance, digits = 2)
names(covTable) <- 1:7
#names(covTable) <- leaves[1:7]
corMat <- cov2cor(covariance)
corTable <- xtable(corMat, digits = 2)
rownames(corTable) <- leaves[1:7]

