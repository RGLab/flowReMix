sqrtMat <- function(X) {
  if(length(X) > 1) {
    return(expm::sqrtm(X))
  } else {
    return(matrix(sqrt(X)))
  }
}

binomDensity <- function(v, prop, N, eta) {
  result <- sum(dbinom(prop, N, expit(eta + v), log = TRUE))
  return(result)
}

logit <- function(p) log(p / (1 - p))

subsetResponseMixtureNested <- function(formula, sub.population = NULL,
                                  N = NULL, id,
                                  data = parent.frame(),
                                  treatment,
                                  weights = NULL,
                                  updateLag = 3,
                                  centerCovariance = TRUE,
                                  rate = 1, nsamp = 100,
                                  maxIter = 30, tol = 1e-03) {
  call <- as.list(match.call())
  if(is.null(call$treatment)) {
    stop("Treatment variable must be specified!")
  }

  #### Getting all relevant variables from call
  # Getting model frame
  dat <- model.frame(formula, data, na.action=na.pass)
  n <- dim(dat)[1]

  # Getting id, waves, weights and treatment variable
  if(typeof(data) == "environment"){
    id <- id
    if(is.null(call$id)) stop("id must be specified!")
    weights <- weights
    if(is.null(call$weights)) weights <- rep(1, n)
    waves <- waves
    if(is.null(call$waves)) {
      warning("No waves variables specified, assuming dataset ordered according to waves!")
      uniqueID <- unique(id)
      idIndex <- lapply(uniqueID,function(x) which(id==x))
      waves <- numeric(length(id))
      for(i in 1:length(idIndex)) {
        waves[idIndex[[i]]] <- 1:length(idIndex[[i]])
      }
    }
    treatment <- treatment
  }
  else{
    if(length(call$id) == 1){
      subj.col <- which(colnames(data) == call$id)
      if(length(subj.col) > 0){
        id <- data[,subj.col]
      }else{
        id <- eval(call$id, envir = parent.frame())
      }
    }else if(is.null(call$id)){
      id <- 1:n
    }

    if(length(call$weights) == 1){
      weights.col <- which(colnames(data) == call$weights)
      if(length(weights.col) > 0){
        weights <- data[, weights.col]
      }else{
        weights <- eval(call$weights, envir=parent.frame())
      }
    }else if(is.null(call$weights)){
      weights <- rep.int(1, n)
    }

    if(length(call$treatment) == 1){
      treatment.col <- which(colnames(data) == call$treatment)
      if(length(treatment.col) > 0){
        treatment <- data[,treatment.col]
      }else{
        treatment <- eval(call$treatment, envir=parent.frame())
      }
    }else if(is.null(call$treatment)){
      stop("Treatment variable must be specified!")
    }
  }

  # getting offset, if no offset is given, then set to zero
  offset <- model.offset(dat)
  if(is.null(offset)){
    off <- rep(0, n)
  }else{
    off <- offset
  }

  # getting N
  if(typeof(data) == "environment"){
    N <- N
  }
  else {
    if(length(call$N) == 1){
      N.col <- which(colnames(data) == call$N)
      if(length(N.col) > 0) {
        N <- data[,N.col]
      }
      else {
        N <- eval(call$N, envir=parent.frame())
      }
    }
    else if(is.null(call$N)){
      stop("N must be specified!")
    }
  }

  # getting Subpopulation
  if(typeof(data) == "environment"){
    sub.population <- sub.population
  }
  else {
    if(length(call$sub.population) == 1){
      s.col <- which(colnames(data) == call$sub.population)
      if(length(s.col) > 0) {
        sub.population <- data[, s.col]
      }
      else {
        sub.population <- eval(call$sub.population, envir=parent.frame())
      }
    }
    else if(is.null(call$sub.population)){
      sub.population <- rep(1, n)
      sub.population <- factor(sub.population)
    }
  }

  # Sub-population must be a factor
  if(!is.factor(sub.population)) stop("Sub-population must be a factor!")

  #################################
  # Creating working dataset
  y <- model.frame(formula, data)
  y <- model.response(y)
  dat$y <- y
  dat$id <- id
  dat$weights <- weights
  dat$N <- N
  dat$off <- off
  dat$sub.population <- sub.population
  dat$treatment <- treatment
  subpopInd <- as.numeric(dat$sub.population)
  uniqueSubpop <- sort(unique(subpopInd))
  dat$subpopInd <- subpopInd
  glmformula <- update.formula(formula, cbind(y, N - y) ~ .  + offset(randomOffset))
  initFormula <- update.formula(formula, cbind(y, N - y) ~ . + (1|id))

  glmFits <- by(dat, dat$sub.population, function(X) lme4::glmer(initFormula, family = "binomial", data = X))
  coefficientList <- lapply(glmFits, function(x) colMeans(coef(x)[[1]]))

  # Estimating covariance structure from marginal fits (step 0)
  nSubjects <- length(unique(dat$id))
  nSubsets <- length(glmFits)
  levelProbs <- rep(0.5, nSubsets)
  estimatedRandomEffects <- do.call("cbind", sapply(glmFits, function(x) lme4::ranef(x)))
  covariance <- cov(estimatedRandomEffects)

  # Computing new posterior porbabilities
  dat$tempTreatment <- dat$treatment
  databyid <- by(dat, dat$id, function(x) x)
  dat$subpopInd <- as.numeric(dat$sub.population)
  sampCoef <- 0.00001
  sampcov <- sampCoef * covariance
  sqrtcov <- sqrtMat(sampcov)
  invcov <- solve(covariance)
  invSampcov <- solve(sampcov)
  covDet <- as.numeric(determinant(covariance, logarithm = TRUE)$modulus)
  sampCovDet <- as.numeric(determinant(sampcov, logarithm = TRUE)$modulus)
  allPosteriors <- matrix(nrow = nSubjects, ncol = maxIter)
  randomSampList <- lapply(1:2, function(x) x)
  posteriors <- matrix(0, nrow = nSubjects, ncol = nSubsets)
  clusterAssignments <- matrix(0.5, nrow = nSubjects, ncol = nSubsets)
  iterCoefMat <- matrix(ncol = length(glmFits), nrow = maxIter + 1)
  MHcoef <- 0.4
  probSamples <- 100
  clusterDensities <- matrix(nrow = 2, ncol = probSamples)
  isingCoefs <- matrix(0, ncol = nSubsets, nrow = nSubsets)
  for(iter in 1:maxIter) {
    # Refitting Model with current means
    dataByPopulation <- data.frame(data.table::rbindlist(databyid))
    dataByPopulation <- by(dataByPopulation, dataByPopulation$sub.population, function(x) x)
    if(iter > 1) {
      glmFits <- lapply(dataByPopulation, function(popdata)
        glm(glmformula, family = "binomial", data = popdata, weights = weights))
      coefficientList <- mapply(function(coef, fit) coef + (coef(fit) - coef)/max(iter - updateLag, 1)^rate,
                                coefficientList, glmFits, SIMPLIFY = FALSE)
    }

    # Updating Prediction
    for(j in 1:length(dataByPopulation)) {
      if(is.factor(dataByPopulation[[j]]$treatment)) {
        dataByPopulation[[j]]$treatment <- factor(0, levels = levels(dataByPopulation[[j]]$tempTreatment))
      } else {
        dataByPopulation[[j]]$treatment <- 0
      }
      offsetBackup <- dataByPopulation[[j]]$randomOffset
      dataByPopulation[[j]]$randomOffset <- 0
      newNullEta <- predict(glmFits[[j]], newdata = dataByPopulation[[j]])
      dataByPopulation[[j]]$treatment <- dataByPopulation[[j]]$tempTreatment
      newAltEta <- predict(glmFits[[j]], newdata = dataByPopulation[[j]])
      nullEta <- ifelse(iter == 1, 0, dataByPopulation[[j]]$nullEta)
      altEta <- ifelse(iter == 1, 0, dataByPopulation[[j]]$altEta)
      dataByPopulation[[j]]$nullEta <- nullEta + (newNullEta - nullEta)/max(iter - updateLag, 1)^rate
      dataByPopulation[[j]]$altEta <- altEta + (newAltEta - altEta)/max(iter - updateLag, 1)^rate
      dataByPopulation[[j]]$randomOffset <- offsetBackup
    }

    databyid <- do.call("rbind", dataByPopulation)
    databyid <- with(databyid, databyid[order(sub.population, id, decreasing = FALSE), ])
    databyid <- by(databyid, databyid$id, function(x) x)

    iterAssignments <- matrix(0, nrow = nSubjects, ncol = nSubsets)
    MHattempts <- 0
    MHsuccess <- 0
    dataLength <- 0
    MHlag <- 5
    condvar <- 1 / diag(invcov)
    assignmentList <- list()
    assignListLength <- 0
    for(i in 1:nSubjects) {
      subjectData <- databyid[[i]]
      popInd <- subjectData$subpopInd
      singlePopInd <- sapply(sort(unique(popInd)), function(x) which(popInd == x)[1])
      N <- subjectData$N
      y <- subjectData$y
      prop <- y/N
      iterPosteriors <- rep(0, nSubsets)

      # Gibbs sampler for cluster assignments
      for(m in 1:nsamp) {
        for(j in 1:nSubsets) {
          nullEta <- subjectData$nullEta[popInd == j]
          altEta <- subjectData$altEta[popInd == j]
          sigmaHat <- sqrt(covariance[j, j])
          empEta <- logit(prop[popInd == j] + 10^-5)
          subsetCount <- y[popInd == j]
          subsetN <- N[popInd == j]

          # Integrating Densities
          for(k in 0:1) {
            if(k == 0) {
              eta <- nullEta
            } else {
              eta <- altEta
            }
            etaResid <- empEta - eta
            muHat <- mean(etaResid)
            vsample <- rnorm(probSamples, mean = muHat, sd = sigmaHat * MHcoef)
            sampNormDens <- dnorm(vsample, mean = muHat, sd = sigmaHat * MHcoef, log = TRUE)
            normDens <- dnorm(vsample, mean = 0, sd = sigmaHat, log = TRUE)
            importanceWeights <- normDens - sampNormDens
            randomEta <- sapply(eta, function(x) x + vsample)
            binomDensity <- colMeans(dbinom(subsetCount, subsetN, t(expit(randomEta)), log = TRUE))
            clusterDensities[k + 1, ] <- binomDensity + importanceWeights
          }

          densityRatio <- rowSums(exp(clusterDensities - max(clusterDensities)))
          if(m >= 2) {
            priorProb <- expit(sum(c(1, clusterAssignments[i, -j]) * isingCoefs[j, ]))
          } else {
            priorProb <- 0.5
          }
          densityRatio <- densityRatio[1] / densityRatio[2] * (1 - priorProb) / priorProb
          pResponder <- 1 / (1 + densityRatio)
          assignment <- rbinom(1, 1, pResponder)
          clusterAssignments[i, j] <- assignment
          if(assignment == 1) iterPosteriors[j] <- iterPosteriors[j] + 1
        }

        if((m %% 5) == 0) {
          assignmentList[[assignListLength + 1]] <- clusterAssignments[i, ]
          assignListLength <- assignListLength + 1
        }
      }

      # Updating global posteriors
      iterPosteriors <- iterPosteriors / nsamp
      print(iterPosteriors)
      currentPost <- posteriors[i, ]
      posteriors[i, ] <- currentPost * (max(iter - updateLag, 1) - 1)/max(iter - updateLag, 1) + iterPosteriors/max(iter - updateLag, 1)
      #posteriors[i, ] <- iterPosteriors

      # MH sampler for random effects
      randomEst <- as.numeric(estimatedRandomEffects[i, ])
      for(m in 1:nsamp) {
        for(j in 1:nSubsets) {
          eta <- nullEta[popInd == j]
          if(clusterAssignments[i, j] == 0) {
            eta <- subjectData$nullEta[popInd == j]
          } else {
            eta <- subjectData$altEta[popInd == j]
          }
          yj <- y[popInd == j]
          Nj <- N[popInd == j]
          MHattempts <- MHattempts + 1
          current <- as.numeric(randomEst[j])
          condmean <- as.numeric(condvar[j] * invcov[j, -j] %*% (randomEst[-j]))
          proposal <- rnorm(1, mean = current, sd = sqrt(covariance[j, j]) * MHcoef)
          newdens <- sum(dbinom(yj, Nj, expit(eta + proposal), log = TRUE))
          newdens <- newdens + dnorm(proposal, mean = condmean, sd = sqrt(condvar[j]))
          olddens <- sum(dbinom(yj, Nj, expit(eta + as.numeric(randomEst[j])), log = TRUE))
          olddens <- olddens + dnorm(current, mean = condmean, sd = sqrt(condvar[j]))
          if(runif(1) < exp(newdens - olddens))  {
            randomEst[j] <- proposal
            MHsuccess <- MHsuccess + 1
          }
        }
        #if(i < 10) print(round(c(it = iter, i = i, m = m, REST = as.numeric(randomEst)), 2))
      }
      # Updating global estimates
      currentRandomEst <- estimatedRandomEffects[i, ]
      estimatedRandomEffects[i, ] <- currentRandomEst +
        (randomEst - currentRandomEst) / max(iter - updateLag, 1)

      # preparing data for glm
      subjectData$randomOffset[1:length(popInd)] <- as.numeric(randomEst[popInd])
      for(j in 1:nSubsets) {
        # Setting treatment according to cluster assignment
        if(clusterAssignments[i, j] == 1) {
          subjectData$treatment[popInd == j] <- subjectData$tempTreatment[popInd == j]
        } else {
          subjectData$treatment[popInd == j] <- 0
        }
      }
      databyid[[i]] <- subjectData
    }

    # Updating Covariance
    if(iter > updateLag) {
      covariance <- cov.wt(estimatedRandomEffects, rep(1, nrow(estimatedRandomEffects)), center = centerCovariance)$cov
      invcov <- solve(covariance)
      #print(round(cov2cor(covariance), 3))
    }

    # Updating ising
    require(glmnet)
    randomizeAssignments <- function(x, prob = 0.5) {
      if(runif(1) < prob) {
        coordinate <- sample.int(length(x), 1)
        x[coordinate] <- ifelse(x[coordinate] == 1, 0, 1)
      }

      return(x)
    }
    assignmentList <- do.call("rbind",assignmentList)
    assignmentList <- t(apply(assignmentList, 1, randomizeAssignments))
    isingfit <- IsingFit::IsingFit(assignmentList, AND = FALSE,
                                   progressbar = FALSE, plot = FALSE)
    newcoefs <- isingfit$weiadj
    for(j in 1:nSubsets) {
      newcoefs[j, ] <- c(newcoefs[j, j], newcoefs[j, -j])
      isingCoefs[j, ] <- isingCoefs[j, ] + (newcoefs[j, ] - isingCoefs[j, ]) / max(1, iter - updateLag)
    }
    print(isingCoefs)

    MHrate <- MHsuccess / MHattempts
    if(MHrate > 0.35) {
      MHcoef <- MHcoef * 1.5
    } else if(MHrate > 0.234) {
      MHcoef <- MHcoef * 1.05
    } else if(MHrate < 0.15) {
      MHcoef <- MHcoef * .5
    } else {
      MHcoef <- MHcoef * .95
    }

    levelProbs <- colMeans(posteriors)
    print(c(iter, levelProbs))
    print(c(iter, round(sapply(coefficientList, function(x) x[2]), 2)))
    print(apply(posteriors, 2, function(x) round(as.numeric(roc(vaccine ~ x)$auc), 3)))
    print(c(MHcoef, MHrate))
  }

  uniqueIDs <- sapply(databyid, function(x) x$id[1])
  result <- list()
  result$levelProbs <- levelProbs
  result$coefficients <- coefficientList
  result$posteriors <- cbind(uniqueIDs, 1 - posteriors)
  result$covariance <- covariance
  result$nullRandomEffects <- cbind(uniqueIDs, estimatedRandomEffects[seq(from = 1,
                                                         to = nrow(estimatedRandomEffects) - 1,
                                                         by = 2), ])
  result$altRandomEffects <- cbind(uniqueIDs, estimatedRandomEffects[seq(from = 2,
                                                                          to = nrow(estimatedRandomEffects),
                                                                          by = 2), ])
  result$isingCov <- isingCoefs
  result$isingfit <- isingfit
  return(result)
}

