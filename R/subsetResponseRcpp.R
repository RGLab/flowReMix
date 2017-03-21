flowReMix_control <- function(updateLag = 5, randomAssignProb = 0.0, nsamp = 20,
                              initMHcoef = 0.4, dataReplicates = NULL, maxDispersion = 10^3,
                              keepEach = 5, centerCovariance = TRUE, intSampSize = 100) {
  object <- list(updateLag = updateLag,
                 randomAssignProb = randomAssignProb,
                 nsamp = nsamp,
                 initMHcoef = initMHcoef,
                 dataReplicates = dataReplicates,
                 maxDispersion = maxDispersion,
                 keepEach = keepEach,
                 centerCovariance = centerCovariance,
                 intSampSize = intSampSize)
  class(object) <- "flowReMix_control"
  return(object)
}

randomizeAssignments <- function(x, prob = 0.5) {
  if(runif(1) < prob) {
    coordinate <- sample.int(length(x), 1)
    x[coordinate] <- ifelse(x[coordinate] == 1, 0, 1)
  }

  return(x)
}

binomDensity <- function(v, prop, N, eta) {
  result <- sum(dbinom(prop, N, expit(eta + v), log = TRUE))
  return(result)
}

logit <- function(p) log(p / (1 - p))

updateCoefs <- function(coef, fit, iter, updateLag, rate, flag) {
  if(flag > 0 & flag <= 2) {
    coef <- c(coef[1], rep(0, length(coef) - 1))
    return(coef)
  }

  if(class(fit)[1] == "cv.glmnet") {
    update <- as.numeric(coef(fit, s = "lambda.min"))
  } else {
    update <- as.numeric(coef(fit))
  }
  notna <- !is.na(update)
  coef[notna] <- coef[notna] + (update[notna] - coef[notna])/max(iter - updateLag, 1)^rate
  return(coef)
}

replicateDataset <- function(data, replicate) {
  data$id <- paste(data$id, "%%%", replicate, sep = "")
  return(data)
}

estimateIntercept <- function(propMat) {
  prop <- pmin(propMat[, 1] + 10^-5, .99)
  estProp <- pmin(propMat[, 2] + 10^-5, .99)
  logit <- log(prop/(1 - prop))
  estLogit <- log(estProp / (1 - estProp))
  return(mean(logit - estLogit))
}

initializeModel <- function(dat, formula) {
  if(is.null(dat)) {
    warning("Some cell-subsets are empty!")
    return("empty!")
  }
  if(all(dat$treatmentvar == 1)) {
    props <- dat$y / dat$N
    mu <- mean(props)
    var <- var(props)
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    probs <- pbeta(props, alpha, beta)
    dat$treatmentvar <- rbinom(nrow(dat), 1, probs)
  }
  X <- model.frame(formula, data = dat)
  X <- model.matrix(formula, X)[, -1]
  fit <- glmnet::cv.glmnet(X, y =  cbind(dat$y, dat$N - dat$y), family = "binomial", weights = dat$weights)
  coef <- glmnet::coef.cv.glmnet(fit, s = "lambda.min")[, 1]
  prop <- dat$y / dat$N
  estProp <- glmnet::predict.cv.glmnet(fit, type = "response", newx = X)
  propMat <- cbind(prop, estProp)
  randomEffects <- as.numeric(by(propMat, dat$id, estimateIntercept))
  randomEffects <- randomEffects[order(unique(dat$id))]
  return(list(fit = fit, coef = coef, randomEffects = randomEffects))
}

flowReMix <- function(formula,
                      cell_type = NULL,
                      subject_id,
                      data = parent.frame(),
                      cluster_variable,
                      cluster_assignment = NULL,
                      weights = NULL,
                      covarianceMethod = c("sparse", "dense", "diagonal"),
                      isingMethod = c("sparse", "dense", "none"),
                      regressionMethod = c("binom", "betabinom", "sparse"),
                      iterations = 10, verbose = TRUE,
                      control = NULL) {
  # Getting control variables -------------------
  if(is.null(control)) {
    control <- flowReMix_control()
  } else if(class(control) != "flowReMix_control") {
    stop("`control' object must be of `flowReMix_control' class!")
  }
  updateLag <- control$updateLag
  randomAssignProb <- control$randomAssignProb
  nsamp <- control$nsamp
  dataReplicates <- control$dataReplicates
  maxDispersion <- control$maxDispersion
  centerCovariance <- control$centerCovariance
  intSampSize <- control$intSampSize
  initMHcoef <- control$initMHcoef
  keepEach <- control$keepEach

  # Checking if inputs are valid --------------------------
  maxIter <- iterations
  rate <- 1
  updateLag <- max(ceiling(updateLag), 1)
  if(updateLag < 2) {
    warning("We recommend using an updateLag of at least 3 to let the algorithm warm up.")
  }

  maxIter <- max(ceiling(maxIter), 1)
  if(maxIter < 4) {
    warning("We recommend running the EM algorithm for at least 4 iterations.")
  }

  nsamp <- ceiling(nsamp)
  if(nsamp < 5) {
    nsamp <- 5
    warning("Number of samples per MH iteration must be 5 or larger!")
  }

  call <- as.list(match.call())
  if(is.null(call$cluster_variable)) {
    stop("Treatment variable must be specified!")
  }

  if(length(isingMethod) > 1) isingModel<- isingModel[1]
  if(!(isingMethod %in% c("sparse", "dense", "none"))) {
    stop("isingMethod must be one of sparse, dense or none")
  }

  if(length(regressionMethod) > 1) regressionMethod <- regressionMethod[1]
  if(!(regressionMethod %in% c(c("binom", "betabinom", "sparse")))) stop("regressionMethod must be one of binom, betabinom or sparse")
  if(regressionMethod == "binom") {
    smallCounts <- FALSE
    betaDispersion <- FALSE
  } else if(regressionMethod == "betabinom") {
    smallCounts <- FALSE
    betaDispersion <- TRUE
  } else {
    smallCounts <- TRUE
    betaDispersion <- FALSE
  }

  if(length(covarianceMethod) > 1) covarianceMethod <- covarianceMethod[1]
  if(!(isingMethod %in% c("sparse", "dense", "none"))) {
    stop("covarianceMethod must be one of sparse, dense or diagonal")
  }

  #### Getting all relevant variables from call --------------------
  # Getting model frame
  dat <- model.frame(formula, data, na.action=na.pass)
  n <- dim(dat)[1]

  # Getting id, waves, weights and treatment variable
  if(typeof(data) == "environment"){
    id <- subject_id
    if(is.null(call$subject_id)) stop("id must be specified!")
    weights <- weights
    if(is.null(call$weights)) weights <- rep(1, n)
    treatmentvar <- cluster_variable
  }
  else{
    if(length(call$subject_id) == 1){
      subj.col <- which(colnames(data) == call$subject_id)
      if(length(subj.col) > 0){
        id <- data[, subj.col]
      }else{
        id <- eval(call$subject_id, envir = parent.frame())
      }
    }else if(is.null(call$subject_id)){
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

    if(length(call$cluster_variable) == 1){
      treatment.col <- which(colnames(data) == call$cluster_variable)
      if(length(treatment.col) > 0){
        treatmentvar <- data[,treatment.col]
      }else{
        treatmentvar <- eval(call$cluster_variable, envir=parent.frame())
      }
    }else if(is.null(call$cluster_variable)){
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

  # getting Subpopulation
  if(typeof(data) == "environment"){
    sub.population <- cell_type
  }
  else {
    if(length(call$cell_type) == 1){
      s.col <- which(colnames(data) == call$cell_type)
      if(length(s.col) > 0) {
        sub.population <- data[, s.col]
      }
      else {
        sub.population <- eval(call$cell_type, envir=parent.frame())
      }
    }
    else if(is.null(call$cell_type)){
      sub.population <- rep(1, n)
      sub.population <- factor(sub.population)
    }
  }

  # Sub-population must be a factor
  if(!is.factor(sub.population)) stop("Sub-population must be a factor!")

  # Creating working dataset ------------------------------
  y <- model.frame(formula, data)
  y <- model.response(y)
  if(!is.matrix(y) | ncol(y) != 2) {
    stop("Response must be a two columns matrix, the first column of which is the
         count of the cell subset and the second is the `parent count - cell count'.")
  }
  N <- rowSums(y)
  y <- y[, 1]
  dat$y <- y
  dat$N <- N
  dat$id <- id
  dat$weights <- weights
  dat$prop <- y / N
  dat$off <- off
  dat$sub.population <- sub.population
  dat$treatmentvar <- treatmentvar
  dat <- dat[, -1]

  # Determining number of data replicates ----------------------
  if(is.null(dataReplicates)) {
    nSubjects <- length(unique(dat$id))
    dataReplicates <- max(ceiling(300 / nSubjects), 2)
  }

  # replicating data set ---------------------------------------
  if(dataReplicates > 1) {
    if(round(dataReplicates) != dataReplicates) warning("dataReplicates rounded to the nearest positive whole number!")
    dataReplicates <- round(dataReplicates)
    dat <- do.call("rbind", lapply(1:dataReplicates, function(x) replicateDataset(dat, x)))
  } else {
    dataReplicates <- 1
  }

  subpopInd <- as.numeric(dat$sub.population)
  uniqueSubpop <- sort(unique(subpopInd))
  dat$subpopInd <- subpopInd

  # Editing formulas ------------------------------
  covariates <- deparse(formula[[3]])
  covariates <- gsub(as.character(call$cluster_variable), "treatmentvar", covariates)
  formula <- update.formula(formula, as.formula(paste(". ~", covariates)))
  glmformula <- update.formula(formula, cbind(y, N - y) ~ .  + offset(randomOffset))
  initFormula <- update.formula(formula, cbind(y, N - y) ~ .)

  # Initializing covariates and random effects------------------
  dataByPopulation <- by(dat, dat$sub.population, function(x) x)
  initialization <- lapply(dataByPopulation, initializeModel, formula = initFormula)
  isEmpty <- sapply(initialization, function(x) x[1] == "empty!")
  if(any(isEmpty)) {
    empty <- names(initialization)[isEmpty]
    newlevels <- levels(sub.population)[!(levels(sub.population) %in% empty)]
    sub.population <- factor(as.character(sub.population), levels = newlevels)
    initialization <- initialization[!isEmpty]
    dat$sub.population <- sub.population
  }
  coefficientList <- lapply(initialization, function(x) x$coef)
  glmFits <- lapply(initialization, function(x) x$fit)
  nSubjects <- length(unique(dat$id))
  nSubsets <- length(glmFits)
  estimatedRandomEffects <- lapply(initialization, function(x) x$randomEffects)
  estimatedRandomEffects <- lapply(estimatedRandomEffects, function(x) {
    if(length(x) < nSubjects) x <- c(x, sample(x, nSubjects - length(x), replace = TRUE))
    return(x)
    })
  estimatedRandomEffects <- do.call("cbind", estimatedRandomEffects)

  # Initializing covariance from diagonal covariance
  levelProbs <- rep(0.5, nSubsets)
  covariance <- diag(apply(estimatedRandomEffects, 2, var))
  invcov <- diag(1 / diag(covariance))
  isingfit <- NULL

  # Setting up preAssignment ----------------------
  if(is.null(cluster_assignment)) {
    preAssignment <- expand.grid(id = unique(dat$id), subset = unique(dat$sub.population))
    preAssignment$assign <- rep(-1, nrow(preAssignment))
  } else {
    preAssignment <- cluster_assignment
    if(ncol(preAssignment) != 3) stop("preAssignment must have 3 columns: id, subset, assignment")
    if(any(!(preAssignment[, 3] %in% c(-1, 0, 1)))) stop("The third column of preAssignment must take values -1, 0 or 1.")
    preAssignment <- data.frame(preAssignment)
    names(preAssignment) <-  c("id", "subset", "assign")
    if(dataReplicates > 1) {
      preAssignment <- do.call("rbind", lapply(1:dataReplicates, function(x) replicateDataset(preAssignment, x)))
    }
    if(nrow(preAssignment) != (nSubsets * nSubjects)) stop("preAssignment must have nSubjects X nSubsets rows.")
    if(any(!(preAssignment[, 1] %in% unique(dat$id)))) stop("The first column of Preassignment must correspond to the id variable.")
  }

  preAssignment <- preAssignment[order(preAssignment$id, preAssignment$subset), ]
  preAssignment <- by(preAssignment, preAssignment$id, function(x) x)

  # More preparations ---------------------------
  dat$tempTreatment <- dat$treatmentvar
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

  flagEquation <- rep(0, nSubsets)

  # Starting analysis -------------------------------
  for(iter in 1:maxIter) {
    # Refitting Model with current random effects/assignments
    dataByPopulation <- data.frame(data.table::rbindlist(databyid))
    dataByPopulation <- by(dataByPopulation, dataByPopulation$sub.population, function(x) x)
    oldM <- M
    if(iter > 1) {
      if(betaDispersion) {
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
            try(glmFits[[j]] <- glm(glmformula, family = "binomial", data = dataByPopulation[[j]], weights = weights))
          } else {
            glmFits[[j]] <- tempfit
            rho <- exp(tempfit$sigma.coefficients)
            M[j] <- max((1 - rho) / rho, maxDispersion)
          }
        }
      } else if(smallCounts) {
        glmFits <- lapply(dataByPopulation, function(popdata) {
          X <- NULL
          try(X <- model.matrix(glmformula, data = popdata)[, - 1], silent = TRUE)
          if(is.null(X)) return(NULL)
          y <- cbind(popdata$N - popdata$y, popdata$y)
          fit <- glmnet::cv.glmnet(X, y, weights = popdata$weights, family = "binomial")
        })
        } else {
        for(j in 1:nSubsets) {
          tempfit <- NULL
          try(tempfit <- glm(glmformula, family = "binomial", data = dataByPopulation[[j]], weights = weights))
          if(!is.null(tempfit)) {
            glmFits[[j]] <- tempfit
          }
        }
      }
      coefficientList <- mapply(updateCoefs, coefficientList, glmFits,
                                iter, updateLag, rate, flagEquation,
                                SIMPLIFY = FALSE)
    }
    M <- oldM + (M - oldM) / max(1, iter - updateLag)

    # Updating Prediction
    for(j in 1:length(dataByPopulation)) {
      coefs <- coefficientList[[j]]
      if(is.factor(dataByPopulation[[j]]$treatmentvar)) {
        dataByPopulation[[j]]$treatmentvar <- factor(0, levels = levels(dataByPopulation[[j]]$tempTreatment))
      } else {
        dataByPopulation[[j]]$treatmentvar <- 0
      }

      modelMat <- NULL
      try(modelMat <- model.matrix(initFormula, data = dataByPopulation[[j]]))
      if(!all(colnames(modelMat) %in% names(coefs))) {
        modelMat <- modelMat[, colnames(modelMat) %in% names(coefs)]
      }
      newNullEta <- as.numeric(modelMat %*% coefs)

      dataByPopulation[[j]]$treatmentvar <- dataByPopulation[[j]]$tempTreatment
      modelMat <- model.matrix(initFormula, data = dataByPopulation[[j]])
      if(!all(colnames(modelMat) %in% names(coefs))) {
        modelMat <- modelMat[, colnames(modelMat) %in% names(coefs)]
      }
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
    # S-step
    for(i in 1:nSubjects) {
      subjectData <- databyid[[i]]
      popInd <- subjectData$subpopInd
      N <- subjectData$N
      y <- subjectData$y
      prop <- y/N
      iterPosteriors <- rep(0, nSubsets)

      # Gibbs sampler for cluster assignments
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
          subjectData$treatmentvar[popInd == j] <- subjectData$tempTreatment[popInd == j]
        } else {
          subjectData$treatmentvar[popInd == j] <- 0
        }
      }
      databyid[[i]] <- subjectData
    }

    # Updating Covariance
    randomList <- do.call("rbind", randomList)
    oldCovariance <- covariance
    if(iter > 1) {
      if(covarianceMethod == "sparse") {
        pdsoftFit <- PDSCE::pdsoft.cv(randomList, init = "dense")
        covariance <- pdsoftFit$sigma
        invcov <- pdsoftFit$omega
      } else if(covarianceMethod == "dense") {
        covariance <- cov.wt(randomList, rep(1, nrow(randomList)), center = centerCovariance)$cov
        invcov <- solve(covariance)
      } else if(covarianceMethod == "diagonal") {
        if(centerCovariance) {
          randomList <- apply(randomList, 2, function(x) x - mean(x))
        }
        covariance <- diag(apply(randomList, 2, function(x) mean(x^2)))
        invcov <- diag(1 / diag(covariance))
      }
    }

    # Updating ising
    if(iter == maxIter) {
      exportAssignment <- assignmentList
      names(exportAssignment) <- names(databyid)
    }
    assignmentList <- do.call("rbind",assignmentList)
    unscrambled <- assignmentList
    if(randomAssignProb > 0 & randomAssignProb <= 1) {
      assignmentList <- t(apply(assignmentList, 1, randomizeAssignments, prob = randomAssignProb))
    }
    assignmentList <- data.frame(assignmentList)
    names(assignmentList) <- names(dataByPopulation)

    if(isingMethod == "sparse" & nSubsets > 2) {
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
      }
    } else if(isingMethod == "dense") {
      for(j in 1:nSubsets) {
        firth <- logistf::logistf(assignmentList[, j] ~ assignmentList[, -j],
                                  datout = FALSE)
        firth <- coef(firth)
        intercept <- firth[1]
        firth[-j] <- firth[-1]
        firth[j] <- intercept
        isingCoefs[j, ] <- isingCoefs[j, ] + (firth - isingCoefs[j, ]) /  max(1, iter - updateLag)
      }
    } else {
      levelProbabilities <- colMeans(assignmentList)
      isingCoefs <- matrix(0, nrow = nSubsets, ncol = nSubsets)
      minprob <- 10^-4
      diag(isingCoefs) <- logit(pmin(pmax(levelProbabilities, minprob), 1 - minprob))
    }
    levelProbs <- colMeans(posteriors)

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
      print(c(iter, levelP = levelProbs))
      print(c(iter, coef = as.numeric(round(sapply(coefficientList, function(x) x[2]), 2))))
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
    names(posteriors) <- c(as.character(call$id), names(dataByPopulation))
  } else {
    realIDs <- gsub("\\%%%.*", "", uniqueIDs)
    post <- by(posteriors, INDICES = realIDs, FUN = colMeans)
    postid <- names(post)
    posteriors <- data.frame(do.call("rbind", post))
    names(posteriors) <- names(dataByPopulation)
    posteriors <- cbind(id = postid, 1 - posteriors)
    names(posteriors)[1] <- as.character(call$id)
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
  result$assignmentList <- exportAssignment
  return(result)
}

