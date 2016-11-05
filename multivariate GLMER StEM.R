#set.seed(502)
set.seed(504)

# data <- rv144
# leaves <- unique(data$population)
# data <- subset(data, population %in% leaves[c(1:7)])
# data <- subset(data, stim != "sebctrl")
# data$treatment <- as.numeric(data$stim == "env")
# data$ptid <- as.numeric(data$ptid)
# data$ptid[data$vaccine == "VACCINE"] <- data$ptid[data$vaccine == "VACCINE"] * 10^4
# data$prop <- data$count / data$parentcount
# data$population <- as.factor(data$population)
# data <- data[order(data$population, data$ptid, data$stim, decreasing = FALSE), ]
#
# mixtureFitList <- by(data, data$population, function(X)
#                      mixedfit <- glmmMixture(count ~ (age + gender + treatment),
#                                              N = parentcount,
#                                              sub.population = population,
#                                              id = ptid,
#                                              treatment = treatment,
#                                              data = X,
#                                              tol = 0.01,
#                                              maxiter = 20,
#                                              nAGQ = 1))
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
logLikelihoods <- matrix(nrow = 2, ncol = nsamp)
effectWeightsMat <- logLikelihoods
posteriors <- numeric(nSubjects)
clusterAssignments <- numeric(nSubjects)
rate <- 1
lastMean <- randomEffects
maxIter <- 70
iterCoefMat <- matrix(ncol = length(mixtureFitList), nrow = maxIter + 1)
posteriorMat <- matrix(ncol = nSubjects, nrow = maxIter)
levelProbsMat <- matrix(ncol = 2, nrow = maxIter)
for(iter in 1:maxIter) {
  nsamp <- 30 + iter
  logLikelihoods <- matrix(nrow = 2, ncol = nsamp)
  zSamp <- matrix(rnorm(nsamp * ncol(randomEffects)), nrow = ncol(randomEffects))
  randomEffectSamp <- sqrtcov %*% zSamp
  accept <- 0
  for(i in 1:nSubjects) {
    subjectData <- databyid[[i]]
    popInd <- subjectData$subpopInd
    N <- subjectData$parentcount
    y <- subjectData$count
    prop <- y/N

    # Performing MH step
    for(k in 1:2) {
      currentSamp <- lastMean[2*i - 2 + k, ]
      sampMean <- 0
      randEst <- estimatedRandomEffects[2*i - 2 + k, ]
      mu <- expit(eta + currentSamp[popInd])
      currentloglik <- sum(dbinom(y, N, mu, log = TRUE)) - 0.5 * t(currentSamp - randEst) %*% invcov %*% (currentSamp - randEst)
      for(j in 1:ncol(randomEffectSamp)) {
        newSamp <- randomEffectSamp[, j] + currentSamp
        mu <- expit(eta + newSamp[popInd])
        newlogLik <- sum(dbinom(y, N, mu, log = TRUE)) - 0.5 * t(newSamp - randEst) %*% invcov %*% (newSamp - randEst)
        if(runif(1) < exp(newlogLik - currentloglik)) {
          currentSamp <- newSamp
          currentloglik <- newlogLik
          accept <- accept + 1
        }
        sampMean <- sampMean + currentSamp / nsamp
      }
      lastMean[2*i - 2 + k, ] <- sampMean
      currentEst <- estimatedRandomEffects[2*i - 2 + k, ]
      estimatedRandomEffects[2*i - 2 + k, ] <- currentEst + (currentSamp - currentEst) / (iter)
    }

    # Computing posteriors and selecting cluster
    for(k in 1:2) {
      #randEst <- estimatedRandomEffects[2*i - 2 + k, ]
      randEst <- lastMean[2*i - 2 + k, ]
      randSamp <- apply(randomEffectSamp, 2, function(x) x + randEst)
      if(iter == 1) {
        if(k == 1) {
          eta <- logit(subjectData$nullMu) - randEst[popInd]
          subjectData$nullEta <- eta
        } else {
          eta <- logit(subjectData$altMu) - randEst[popInd]
          subjectData$altEta <- eta
        }
      } else {
        if(k == 1) {
          eta <- subjectData$nullEta
        } else {
          eta <- subjectData$altEta
        }
      }

      mu <- expit(apply(randSamp, 2, function(x) eta + x[popInd]))
      binomLlik <- apply(mu, 2, function(x) sum(dbinom(y, N, x, log = TRUE)))
      normWeights <- apply(randSamp, 2, function(x) - 0.5 * t(x) %*% invcov  %*% x - 0.5 * covDet)
      importanceWeights <- apply(randSamp, 2, function(x) - 0.5 * t(x - randEst) %*% invSampcov  %*% (x - randEst) - 0.5 * sampCovDet)
      logLikelihoods[k, ] <- binomLlik + normWeights - importanceWeights + log(levelProbs[k])
      randomSampList[[k]] <- randSamp
    }

    posterior <- rowSums(exp(logLikelihoods - max(logLikelihoods)))
    posterior <-  1/(1 + posterior[1] / posterior[2])
    if(is.nan(posterior)) posterior <- 0
    #print(c(i, subjectData$vaccine[1] == "VACCINE", posterior))
    posteriorMat[iter, i] <- posterior
    if(iter == 1) {
      posteriors[i] <- posterior
    } else if(iter > 1) {
      posteriors[i] <- posteriors[i] + (posterior - posteriors[i]) / (iter)
    }

    # Sampling cluster assignment (Step 1b)
    cluster <- 1 + rbinom(1, 1, posterior)
    clusterAssignments[i] <- cluster
    subjectData$tempTreatment <- subjectData$treatment * (cluster == 2)

    subjectData$randomOffset <- lastMean[2*i - 2 + cluster, ]
    databyid[[i]] <- subjectData
  }

  # Refitting Model with current means
  dataForGlm <- do.call("rbind", databyid)
  dataByPopulation <- by(dataForGlm, dataForGlm$population, function(x) x)
  glmFits <- lapply(dataByPopulation, function(popdata)
                glm(cbind(count, parentcount - count) ~  age + gender + tempTreatment + offset(randomOffset),
                    family = "binomial", data = popdata))
  coefficientList <- mapply(function(coef, fit) coef + (coef(fit) - coef)/(iter),
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
    dataByPopulation[[j]]$nullEta <- newNullEta
    dataByPopulation[[j]]$altEta <- newAltEta
  }

  databyid <- do.call("rbind", dataByPopulation)
  databyid <- with(databyid, databyid[order(population, ptid, stim, decreasing = FALSE), ])
  databyid <- by(databyid, databyid$ptid, function(x) x)

  # Estimating Covariance
  weights <- as.vector(sapply(posteriors, function(x) c(1 - x, x)))
  covariance <- cov.wt(estimatedRandomEffects, weights, center = FALSE)$cov
  invcov <- solve(covariance)
  covDet <- as.numeric(determinant(covariance, logarithm = TRUE)$modulus)

  acceptRate <- accept / (nsamp * nSubjects)
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

  levelProbs[2] <- mean(posteriorMat[iter, ])
  levelProbs[1] <- 1 - levelProbs[2]
  levelProbsMat[iter, ] <- levelProbs

  # Some Diagnostics/Outputs
  iterCoef <- sapply(glmFits, function(x) coef(x)[5])
  iterCoefMat[iter + 1, ] <- iterCoef
  currentCoef <- sapply(coefficientList, function(x) x[[5]])
  initCoef <- sapply(mixtureFitList, function(x) x$beta[[5]])
  iterCoefMat[1, ] <- initCoef
  require(pROC)
  rocfit <- roc(vaccine ~ posteriors)
  print(c(iter, sampCoef, accept / (nsamp * nSubjects)))
  print(colMeans(levelProbsMat[1:iter, , drop = FALSE]))
  print(round(cbind(iterCoef, currentCoef, initCoef), 3))
  #print(round(cov2cor(covariance), 3))
  print(plot(rocfit, main = round(rocfit$auc, 3)))
}


posteriors <- colMeans(posteriorMat[1:maxIter, ])
levelProbs <- colMeans(levelProbsMat[1:maxIter, ])


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
  print(ggplot() +
    geom_point(aes(x = ctrlprop, y = stimprop, col = posteriors, shape = !vaccine, size = !vaccine), fill = "white") +
    theme_bw() + geom_abline(intercept = 0, slope = 1) +
    scale_colour_gradientn(colours=rainbow(4)))
}

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
abline(h = c(.9, .95), v = c(.05, 0.1), col = "grey")
