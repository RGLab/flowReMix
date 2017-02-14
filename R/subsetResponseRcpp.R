randomizeAssignments <- function(x, prob = 0.5) {
  if(runif(1) < prob) {
    coordinate <- sample.int(length(x), 1)
    x[coordinate] <- ifelse(x[coordinate] == 1, 0, 1)
  }

  return(x)
}

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

updateCoefs <- function(coef, fit, iter, updateLag, rate) {
  update <- coef(fit)
  notna <- !is.na(update)
  coef[notna] <- coef[notna] + (update[notna] - coef[notna])/max(iter - updateLag, 1)^rate
  return(coef)
}

replicateDataset <- function(data, replicate) {
  data$id <- paste(data$id, "%%%", replicate, sep = "")
  return(data)
}

subsetResponseMixtureRcpp <- function(formula, sub.population = NULL,
                                  N = NULL, id,
                                  data = parent.frame(),
                                  treatment,
                                  preAssignment = NULL,
                                  weights = NULL,
                                  centerCovariance = TRUE,
                                  covarianceMethod = c("dense" , "sparse"),
                                  sparseGraph = FALSE,
                                  betaDispersion = TRUE,
                                  randomAssignProb = 0.0,
                                  updateLag = 3, rate = 1, nsamp = 100,
                                  initMHcoef = 0.5,
                                  maxIter = 8, verbose = TRUE,
                                  dataReplicates = 1) {
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
  dat$prop <- y / N
  dat$off <- off
  dat$sub.population <- sub.population
  dat$treatment <- treatment
  # replicating
  if(dataReplicates > 1) {
    dat <- do.call("rbind", lapply(1:dataReplicates, function(x) replicateDataset(dat, x)))
  }

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

  # Sometime the initalization will yield a vector of `zero` random effects
  invalid <- apply(estimatedRandomEffects, 2, function(x) all(x == 0))
  if(any(invalid)) {
    estimatedRandomEffects[, invalid] <- rnorm(sum(invalid) * nrow(estimatedRandomEffects), sd = sd(unlist(estimatedRandomEffects)))
  }
  covariance <- PDSCE::pdsoft.cv(estimatedRandomEffects, init = "diag", nsplits = 10)
  invcov <- covariance$omega
  covariance <- covariance$sigma
  isingfit <- NULL

  # Setting up preAssignment
  if(is.null(preAssignment)) {
    preAssignment <- expand.grid(id = unique(dat$id), subset = unique(dat$sub.population))
    preAssignment$assign <- rep(-1, nrow(preAssignment))
  } else {
    if(ncol(preAssignment) != 3) stop("preAssignment must have 3 columns: id, subset, assignment")
    if(nrow(preAssignment) != (nSubsets * nSubjects)) stop("preAssignment must have nSubjects X nSubsets rows.")
    if(any(!(preAssignment[, 3] %in% c(-1, 0, 1)))) stop("The third column of preAssignment must take values -1, 0 or 1.")
    if(any(!(preAssignment[, 1] %in% unique(dat$id)))) stop("The first column of Preassignment must correspond to the id variable.")
    preAssignment <- data.frame(preAssignment)
    names(preAssignment) <-  c("id", "subset", "assign")
  }
  preAssignment <- preAssignment[order(preAssignment$id, preAssignment$subset), ]
  preAssignment <- by(preAssignment, preAssignment$id, function(x) x)

  # More set up....
  dat$tempTreatment <- dat$treatment
  databyid <- by(dat, dat$id, function(x) x)
  dat$subpopInd <- as.numeric(dat$sub.population)
  posteriors <- matrix(0, nrow = nSubjects, ncol = nSubsets)
  clusterAssignments <- matrix(0.5, nrow = nSubjects, ncol = nSubsets)
  iterCoefMat <- matrix(ncol = length(glmFits), nrow = maxIter + 1)
  if(length(initMHcoef) == nSubsets) {
    MHcoef <- initMHcoef
  } else if(is.null(initMHcoef)) {
    MHcoef <- rep(0.5, nSubsets)
  } else {
    MHcoef <- rep(initMHcoef[1], nSubsets)
  }
  probSamples <- 100
  clusterDensities <- matrix(nrow = 2, ncol = probSamples)
  isingCoefs <- matrix(0, ncol = nSubsets, nrow = nSubsets)
  if(betaDispersion) {
    M <- rep(10^4, nSubsets)
  } else {
    M <- rep(10^8, nSubsets)
  }
  for(iter in 1:maxIter) {
    # Refitting Model with current means
    dataByPopulation <- data.frame(data.table::rbindlist(databyid))
    dataByPopulation <- by(dataByPopulation, dataByPopulation$sub.population, function(x) x)
    oldM <- M
    if(iter > 1) {
      if(betaDispersion & iter > 1) {
        for(j in 1:nSubsets) {
          popdata <- dataByPopulation[[j]]
          if(class(glmFits[[j]])[[1]] == "gamlss") {
            startFrom <- glmFits[[j]]
          } else {
            startFrom <- NULL
          }
          tempfit <- NULL
          try(tempfit <- gamlss::gamlss(formula = glmformula, weights = weights,
          family = gamlss.dist::BB, data = popdata, start.from = startFrom))
          if(is.null(tempfit)) {
            try(glmFits[[j]] <- glm(glmformula, family = "binomial",
                                data = dataByPopulation[[j]], weights = weights))
          } else {
            glmFits[[j]] <- tempfit
            rho <- exp(tempfit$sigma.coefficients)
            M[j] <- max((1 - rho) / rho, 10^3)
          }
        }
      } else {
        glmFits <- lapply(dataByPopulation, function(popdata)
          glm(glmformula, family = "binomial", data = popdata, weights = weights))
      }
      coefficientList <- mapply(updateCoefs, coefficientList, glmFits,
                                iter, updateLag, rate, SIMPLIFY = FALSE)
    }
    M <- oldM + (M - oldM) / max(1, iter - updateLag)

    # Updating Prediction
    for(j in 1:length(dataByPopulation)) {
      coefs <- coefficientList[[j]]
      if(is.factor(dataByPopulation[[j]]$treatment)) {
        dataByPopulation[[j]]$treatment <- factor(0, levels = levels(dataByPopulation[[j]]$tempTreatment))
      } else {
        dataByPopulation[[j]]$treatment <- 0
      }
      modelMat <- NULL
      try(modelMat <- model.matrix(formula, data = dataByPopulation[[j]]))
      newNullEta <- as.numeric(modelMat %*% coefs)

      dataByPopulation[[j]]$treatment <- dataByPopulation[[j]]$tempTreatment
      modelMat <- model.matrix(formula, data = dataByPopulation[[j]])
      newAltEta <- as.numeric(modelMat %*% coefs)

      if(iter > 1) {
        nullEta <- dataByPopulation[[j]]$nullEta
        altEta <- dataByPopulation[[j]]$altEta
      } else {
        nullEta <- 0
        altEta <- 0
      }
      dataByPopulation[[j]]$nullEta <- nullEta + (newNullEta - nullEta)/max(iter - updateLag, 1)^rate
      dataByPopulation[[j]]$altEta <- altEta + (newAltEta - altEta)/max(iter - updateLag, 1)^rate
    }

    databyid <- do.call("rbind", dataByPopulation)
    databyid <- with(databyid, databyid[order(sub.population, id, decreasing = FALSE), ])
    databyid <- by(databyid, databyid$id, function(x) x)

    iterAssignments <- matrix(0, nrow = nSubjects, ncol = nSubsets)
    dataLength <- 0
    MHlag <- 5
    condvar <- 1 / diag(invcov)
    assignmentList <- list()
    randomList <- list()
    assignListLength <- 0
    MHattempts <- rep(0, nSubsets)
    MHsuccess <- rep(0, nSubsets)
    for(i in 1:nSubjects) {
      #print(i)
      subjectData <- databyid[[i]]
      popInd <- subjectData$subpopInd
      N <- subjectData$N
      y <- subjectData$y
      prop <- y/N
      iterPosteriors <- rep(0, nSubsets)
      keepEach <- 1
      intSampSize <- 100

      # Gibbs sampler for cluster assignments
      #set.seed(iter  + i * 10^4)
      unifVec <- runif(nsamp * nSubsets)
      normVec <- rnorm(intSampSize)
      assignmentMat <- subsetAssignGibbs(y, prop, N, isingCoefs,
                                         subjectData$nullEta, subjectData$altEta,
                                         covariance, nsamp, nSubsets, keepEach, intSampSize,
                                         MHcoef,
                                         as.integer(popInd),
                                         unifVec, normVec,
                                         M, betaDispersion,
                                         as.integer(preAssignment[[i]]$assign))

      # Updating global posteriors
      iterPosteriors <- colMeans(assignmentMat)
      posteriors[i, ] <- posteriors[i, ] + (iterPosteriors - posteriors[i, ])/max(iter - updateLag, 1)
      clusterAssignments[i, ] <- assignmentMat[nrow(assignmentMat), ]
      assignmentMat <- assignmentMat[seq(from = 5, to = nrow(assignmentMat), by = 5), ]
      assignmentList[[i]] <- assignmentMat

      # MH sampler for random effects
      unifVec <- runif(nsamp * nSubsets)
      eta <- subjectData$nullEta
      responderSubset <- popInd %in% which(clusterAssignments[i, ] == 1)
      eta[responderSubset] <- subjectData$altEta[responderSubset]
      assignment <- as.vector(clusterAssignments[i, ])
      randomEst <- as.numeric(estimatedRandomEffects[i, ])
      randomMat <- randomEffectCoordinateMH(y, N,
                                            as.integer(i), nsamp, nSubsets,
                                            MHcoef,
                                            as.vector(assignment),
                                            as.integer(popInd),
                                            as.numeric(eta),
                                            randomEst,
                                            as.numeric(condvar), covariance, invcov,
                                            MHattempts, MHsuccess,
                                            unifVec,
                                            M, betaDispersion,
                                            keepEach)
      randomEst <- randomMat[nrow(randomMat), ]
      randomList[[i]] <- randomMat


      # Updating global estimates
      currentRandomEst <- estimatedRandomEffects[i, ]
      estimatedRandomEffects[i, ] <- currentRandomEst + (colMeans(randomMat) - currentRandomEst) / max(iter - updateLag, 1)

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
    randomList <- do.call("rbind", randomList)
    oldCovariance <- covariance
    if(iter > 1) {
      if(covarianceMethod[1] == "sparse") {
        pdsoftFit <- PDSCE::pdsoft.cv(randomList, init = "dense")
        covariance <- pdsoftFit$sigma
        invcov <- pdsoftFit$omega
      } else {
        covariance <- cov.wt(randomList, rep(1, nrow(randomList)), center = centerCovariance)$cov
        invcov <- solve(covariance)
      }
    }

    # Updating ising
    assignmentList <- do.call("rbind",assignmentList)
    unscrambled <- assignmentList
    if(randomAssignProb > 0 & randomAssignProb <= 1) {
      assignmentList <- t(apply(assignmentList, 1, randomizeAssignments, prob = randomAssignProb))
    }
    assignmentList <- data.frame(assignmentList)
    names(assignmentList) <- unique(sub.population)

    if(sparseGraph == TRUE & nSubsets > 2) {
      tempfit <- NULL
      try(tempfit <- IsingFit::IsingFit(assignmentList, AND = FALSE,
                                     progressbar = FALSE, plot = FALSE))
      if(!is.null(tempfit)) {
        isingfit <- tempfit
        isingtemp <- isingfit$weiadj
        diag(isingtemp) <- isingfit$thresholds
        isingtemp[isingtemp == Inf] <- 0
        isingtemp[isingtemp == -Inf] <- 0
        isingCoefs <- isingCoefs + (isingtemp - isingCoefs) / max(1, iter - updateLag)
        levelProbs <- colMeans(posteriors)
      } else {
        print("hello")
      }
    } else {
      for(j in 1:nSubsets) {
        firth <- logistf::logistf(assignmentList[, j] ~ assignmentList[, -j],
                                  datout = FALSE)
        firth <- coef(firth)
        intercept <- firth[1]
        firth[-j] <- firth[-1]
        firth[j] <- intercept
        isingCoefs[j, ] <- isingCoefs[j, ] + (firth - isingCoefs[j, ]) /  max(1, iter - updateLag)
      }
    }

    # Updating MH coefficient
    ratevec <- numeric(nSubsets)
    for(j in 1:nSubsets) {
      MHrate <- MHsuccess[j] / MHattempts[j]
      ratevec[j] <- MHrate
      if(MHrate > .285) {
        MHcoef[j] <- MHcoef[j] * 1.5
      } else if(MHrate > 0.234) {
        MHcoef[j] <- MHcoef[j] * 1.1
      } else if(MHrate < .185) {
        MHcoef[j] <- MHcoef[j] * .5
      } else {
        MHcoef[j] <- MHcoef[j] * .90
      }
    }


    if(verbose) {
      #print(isingCoefs)
      print(c(iter, levelP = levelProbs))
      print(c(iter, coef = as.numeric(round(sapply(coefficientList, function(x) x[2]), 2))))
      #try(print(c(AUC = apply(posteriors, 2, function(x) round(as.numeric(roc(!vaccine ~ x)$auc), 3)))))
      print(c(M = M))
      print(round(c(MH = MHcoef), 3))
      print(round(ratevec, 3))
    }
  }

  # Processing posteriors
  posteriors <- data.frame(posteriors)
  uniqueIDs <- sapply(databyid, function(x) x$id[1])
  if(dataReplicates <= 1) {
    posteriors <- cbind(id = uniqueIDs, 1 - posteriors)
    names(posteriors) <- unique(dat$sub.population)
  } else {
    realIDs <- gsub("\\%%%.*", "", uniqueIDs)
    post <- by(posteriors, INDICES = realIDs, FUN = colMeans)
    postid <- names(post)
    posteriors <- data.frame(do.call("rbind", post))
    names(posteriors) <- unique(dat$sub.population)
    posteriors <- cbind(id = postid, 1 - posteriors)
  }

  # Preparing output
  result <- list()
  result$posteriors <- posteriors
  result$levelProbs <- levelProbs
  result$coefficients <- coefficientList
  result$covariance <- covariance
  result$randomEffects <- estimatedRandomEffects
  result$isingCov <- isingCoefs
  result$isingfit <- isingfit
  result$dispersion <- M
  return(result)
}

