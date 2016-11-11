flowRegressionMixture <- function(formula, sub.population = NULL,
                                  N = NULL, id,
                                  data = parent.frame(),
                                  treatment,
                                  weights = NULL,
                                  rate = 1,
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
  subpopInd <- as.numeric(dat$sub.population)
  uniqueSubpop <- sort(unique(subpopInd))
  dat$subpopInd <- subpopInd
  glmformula <- update.formula(formula, cbind(y, N - y) ~ .  + offset(randomOffset))

  mixtureFitList <- by(dat, dat$sub.population, function(X)
    glmmMixture(formula,
                N = N,
                id = id,
                treatment = treatment,
                data = X,
                tol = 0.01,
                maxiter = 1,
                nAGQ = 1))

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

  # Computing new posterior porbabilities
  dat$tempTreatment <- dat$treatment
  muMat <- lapply(mixtureFitList, function(x) x$mu[, 4:5])
  muMat <- do.call("rbind", muMat)
  dat$nullMu <- muMat$nullMu
  dat$altMu <- muMat$altMu
  databyid <- by(dat, dat$id, function(x) x)
  dat$subpopInd <- as.numeric(data$population)
  nSubjects <- length(unique(dat$id))
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
      N <- subjectData$N
      y <- subjectData$y
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
        mu <- expit(eta+randSamp[popInd,])
        binomLlik <- colSums(dbinom(y,N,mu,log=TRUE))
        normWeights <- -0.5*diag(crossprod(crossprod(invcov,randSamp),randSamp))-0.5*covDet #faster than above
        importanceWeights <- -0.5*diag(crossprod(crossprod(invSampcov,randSamp-randEst),randSamp-randEst))-0.5*sampCovDet #faster than above
        logLikelihoods[k, ] <- binomLlik + normWeights - importanceWeights + log(levelProbs[k])
        randomSampList[[k]] <- randSamp
      }

      posterior <- rowSums(exp(logLikelihoods - max(logLikelihoods)))
      posterior <-  1/(1 + posterior[1] / posterior[2])
      if(is.nan(posterior)) posterior <- 0
      if(iter == 1) {
        posteriors[i] <- posterior
      } else if(iter > 1) {
        posteriors[i] <- posteriors[i] + (posterior - posteriors[i]) / (iter)^rate
      }

      # Sampling cluster assignment
      cluster <- 1 + rbinom(1, 1, posteriors[i])
      clusterAssignments[i] <- cluster
      subjectData$treatment <- subjectData$tempTreatment * (cluster == 2)

      # Performing MH step
      unifs = runif(nsamp)
      flowReMix:::MH(lastMean, estimatedRandomEffects, y, N, randomEffectSamp, i, popInd, invcov, accept, iter,rate, unifs, nullEta, altEta);

      subjectData$randomOffset <- lastMean[2*i - 2 + cluster, ]
      databyid[[i]] <- subjectData
    }

    # Refitting Model with current means
    dataForGlm <- data.frame(data.table::rbindlist(databyid))
    dataByPopulation <- by(dataForGlm, dataForGlm$sub.population, function(x) x)
    glmFits <- lapply(dataByPopulation, function(popdata)
      glm(glmformula, family = "binomial", data = popdata, weights = weights))
    coefficientList <- mapply(function(coef, fit) coef + (coef(fit) - coef)/(iter)^rate,
                              coefficientList, glmFits, SIMPLIFY = FALSE)
    # Updating Prediction
    for(j in 1:length(dataByPopulation)) {
      dataByPopulation[[j]]$treatment <- 0
      dataByPopulation[[j]]$randomOffset <- 0
      newNullEta <- predict(glmFits[[j]], newdata = dataByPopulation[[j]])
      dataByPopulation[[j]]$treatment <- dataByPopulation[[j]]$tempTreatment
      newAltEta <- predict(glmFits[[j]], newdata = dataByPopulation[[j]])
      nullEta <- dataByPopulation[[j]]$nullEta
      altEta <- dataByPopulation[[j]]$altEta
      dataByPopulation[[j]]$nullEta <- nullEta + (newNullEta - nullEta)/(iter)^rate
      dataByPopulation[[j]]$altEta <- altEta + (newAltEta - altEta)/(iter)^rate
    }

    databyid <- do.call("rbind", dataByPopulation)
    databyid <- with(databyid, databyid[order(sub.population, id, decreasing = FALSE), ])
    databyid <- by(databyid, databyid$id, function(x) x)

    acceptRate <- accept / (2 * nsamp * nSubjects)
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
    require(pROC)
    rocfit <- roc(vaccines ~ posteriors)
    print(plot(rocfit, main = round(rocfit$auc, 3)))
    print(cov2cor(covariance))
  }

  uniqueIDs <- sapply(databyid, function(x) x$id[1])
  result <- list()
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
