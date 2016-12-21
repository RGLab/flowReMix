sqrtMat <- function(X) {
  if(length(X) > 1) {
    return(expm::sqrtm(X))
  } else {
    return(matrix(sqrt(X)))
  }
}

subsetResponseMixture <- function(formula, sub.population = NULL,
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
  posteriors <- matrix(0.5, nrow = nSubjects, ncol = nSubsets)
  clusterAssignments <- matrix(nrow = nSubjects, ncol = nSubsets)
  iterCoefMat <- matrix(ncol = length(glmFits), nrow = maxIter + 1)
  accept <- 0
  for(iter in 1:maxIter) {
    # Refitting Model with current means
    dataForGlm <- data.frame(data.table::rbindlist(databyid))
    dataByPopulation <- by(dataForGlm, dataForGlm$sub.population, function(x) x)
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
      dataByPopulation[[j]]$randomOffset <- 0
      newNullEta <- predict(glmFits[[j]], newdata = dataByPopulation[[j]])
      dataByPopulation[[j]]$treatment <- dataByPopulation[[j]]$tempTreatment
      newAltEta <- predict(glmFits[[j]], newdata = dataByPopulation[[j]])
      nullEta <- ifelse(iter == 1, 0, dataByPopulation[[j]]$nullEta)
      altEta <- ifelse(iter == 1, 0, dataByPopulation[[j]]$altEta)
      dataByPopulation[[j]]$nullEta <- nullEta + (newNullEta - nullEta)/max(iter - updateLag, 1)^rate
      dataByPopulation[[j]]$altEta <- altEta + (newAltEta - altEta)/max(iter - updateLag, 1)^rate
    }

    databyid <- do.call("rbind", dataByPopulation)
    databyid <- with(databyid, databyid[order(sub.population, id, decreasing = FALSE), ])
    databyid <- by(databyid, databyid$id, function(x) x)

    iterAssignments <- matrix(0, nrow = nSubjects, ncol = nSubsets)
    for(i in 1:nSubjects) {
      # Computing Posterior Probabilities (Step 1a)
      subjectData <- databyid[[i]]
      popInd <- subjectData$subpopInd
      singlePopInd <- sapply(sort(unique(popInd)), function(x) which(popInd == x)[1])
      N <- subjectData$N
      y <- subjectData$y
      prop <- y/N

      # Initializing the random effects
      if(iter == 1) {
        subjectData$randomOffset <- t(estimatedRandomEffects[i, popInd])
        subjectData$nullOffset <- subjectData$randomOffset
        subjectData$altOffset <- subjectData$randomOffset
      }

      # Sampling nsamp times per subject
      for(m in 1:nsamp) {
        # Iterating over populations
        for(j in 1:nSubsets) {
          # Sampling assignment
          # The probability for assignments is proportional to the densities of the subset
          # At the current 'sampled' step. For now, assume independence between assignments.
          popCount <- y[popInd == j]
          popN <- N[popInd == j]
          nullEta <- subjectData$nullEta[popInd == j]
          altEta <- subjectData$altEta[popInd == j]
          nullEffect <- subjectData$nullOffset[popInd == j]
          altEffect <- subjectData$altOffset[popInd == j]
          altProb <- levelProbs[j]

          nulldens <- sum(dbinom(popCount, popN, expit(nullEta + nullEffect), log = TRUE)) + log(1 - altProb)
          altdens <- sum(dbinom(popCount, popN, expit(altEta + altEffect), log = TRUE)) + log(altProb)
          posteriorProb <- 1/(1 + exp(nulldens - altdens))
          popCluster <- rbinom(1, 1, posteriorProb)
          clusterAssignments[i, j] <- popCluster
          sampledEta <- (popCluster == 1) * altEta + (popCluster == 0) * nullEta
          subjectData$sampledEta[popInd == j] <- (popCluster == 1) * altEta + (popCluster == 0) * nullEta
          subjectData$randomOffset[popInd == j] <- (popCluster == 1) * altEffect + (popCluster == 0) * nullEffect

          # Sample random effect, but only if all cluster assignments are availble.
          randomEffects <- subjectData$randomOffset[singlePopInd]
          if(iter > 1 | m > 1) {
            condvar <- 1 / invcov[j, j]
            #condmean <- -condvar * invcov[j, -j] %*% (randomEffects[-j])
            #proposal <- rnorm(1, mean = condmean, sd = sqrt(condvar))
            proposal <- rnorm(1, mean = randomEffects[j], sd = sqrt(condvar))

            # Because of the gibbs form of the proposal, the MHratio is just the likelihood ratio.
            olddens <- sum(dbinom(popCount, popN, expit(sampledEta + randomEffects[j]), log = TRUE))
            newdens <- sum(dbinom(popCount, popN, expit(sampledEta + proposal), log = TRUE))
            MHratio <- exp(newdens - olddens)
            if(runif(1) < MHratio) {
              subjectData$randomOffset[popInd == j] <- proposal
              if(popCluster == 1) {
                subjectData$altOffset[popInd == j] <- proposal
              } else {
                subjectData$nullOffset[popInd == j] <- proposal
              }
            }
          }
        }
        iterAssignments[i, ] <- iterAssignments[i, ] + clusterAssignments[i, ] / nsamp
      }

      # Setting treatment according to cluster assignment
      for(j in 1:nSubsets) {
        if(clusterAssignments[i, j] == 1) {
          subjectData$treatment[popInd == j] <- subjectData$tempTreatment[popInd == j]
        } else {
          subjectData$treatment[popInd == j] <- 0
        }
      }

      # Updating estimated random effects and posteriors
      currentRandomEst <- estimatedRandomEffects[i, ]
      estimatedRandomEffects[i, ] <- currentRandomEst +
        (subjectData$randomOffset[singlePopInd] - currentRandomEst) / max(iter - updateLag, 1)^rate

      databyid[[i]] <- subjectData
    }

    if(iter > updateLag) {
      covariance <- cov.wt(estimatedRandomEffects, rep(1, nrow(estimatedRandomEffects)), center = centerCovariance)$cov
      invcov <- solve(covariance)
    }
    print(round(covariance, 3))

    # Updating posteriors
    if(iter > updateLag) {
      posteriors <- posteriors * (iter - updateLag - 1) / (iter - updateLag) + iterAssignments / (iter - updateLag)
    } else {
      posteriors <- iterAssignments
    }

    levelProbs <- colMeans(posteriors)
    print(c(iter, levelProbs))
    print(c(iter, round(sapply(coefficientList, function(x) x[2]), 2)))
  }

  uniqueIDs <- sapply(databyid, function(x) x$id[1])
  result <- list()
  result$levelProbs <- levelProbs
  result$coefficients <- coefficientList
  result$posteriors <- cbind(uniqueIDs, 1 - posteriors, posteriors)
  result$covariance <- covariance
  result$nullRandomEffects <- cbind(uniqueIDs, estimatedRandomEffects[seq(from = 1,
                                                         to = nrow(estimatedRandomEffects) - 1,
                                                         by = 2), ])
  result$altRandomEffects <- cbind(uniqueIDs, estimatedRandomEffects[seq(from = 2,
                                                                          to = nrow(estimatedRandomEffects),
                                                                          by = 2), ])
  return(result)
}

