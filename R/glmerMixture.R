# A function for monte-carlo computation of mixture densities
simulateGmmDensity <- function(x, randomsd, nsamp, levelProbs) {
  sampmean <- x$randomMean[[1]]
  normsamp <- rnorm(nsamp, sampmean, randomsd / 2)
  etasamp <- sapply(normsamp, function(v) v + x$eta)
  musamp <- expit(etasamp)
  normdens <- dnorm(normsamp, sd = randomsd, log = TRUE) -
    dnorm(normsamp, mean = sampmean, sd = randomsd / 2)
  binomdens <- apply(musamp, 2, function(v) sum(dbinom(x$count, x$N, v, log = TRUE)))
  levelProb <- log(levelProbs[x$subpopInd[1], x$treatmentLevel[1] + 1])
  return(binomdens + normdens + levelProb)
}

# A function for computation of weights for mixture-GLMER model
computeGmmWeights <- function(dat, randomsd, nsamp = 20, levelProbs) {
  if(all(dat$treatmentLevel == 0)) return(rep(1, nrow(dat)))

  nrows <- nrow(dat)
  eta <- dat$fixedEta
  y <- dat$count
  N <- dat$N
  dat$randomMean <- round(dat$randomMean, 10)
  densities <- by(dat, dat$treatmentLevel, simulateGmmDensity, randomsd, nsamp, levelProbs)
  densities <- do.call("rbind", densities)
  densities <- exp(densities - max(densities))
  densities <- rowSums(densities)
  weights <- numeric(nrow(dat))
  altprob <- 1/(1 + densities[1] / densities[2])
  weights[dat$treatmentLevel == 2] <- altprob
  weights[dat$treatmentLevel == 1] <- 1 - altprob
  return(weights)
}

# Estimate a composite mixture GEE
glmmMixture <- function(formula, sub.population = NULL,
                                 N = NULL, id, waves,
                                 data = parent.frame(),
                                 treatment,
                                 treatment.levels = 2,
                                 compute.offset = TRUE,
                                 weights = NULL,
                                 init.beta = NULL,
                                 init.rho = 0,
                                 rho.fixed = FALSE,
                                 shrinkage.folds = 5,
                                 maxiter = 100, tol = 1e-03,
                                 nAGQ = 1,
                                 computeT = FALSE,
                                 shrink.sandwich = TRUE) {
  #### Some checks
  call <- as.list(match.call())
  if(is.null(call$treatment)) {
    stop("Treatment variable must be specified!")
  }
  if(treatment.levels != round(treatment.levels) | treatment.levels < 1) {
    stop("treatment.levels must be a positive integer > 1.")
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

    if(length(call$waves) == 1){
      waves.col <- which(colnames(data) == call$waves)
      if(length(waves.col) > 0){
        waves <- data[, waves.col]
      }else{
        waves <- eval(call$waves, envir=parent.frame())
      }
    } else if(is.null(call$waves)){
      warning("No waves variables specified, assuming dataset ordered according to waves!")
      uniqueID <- unique(id)
      idIndex <- lapply(uniqueID,function(x) which(id == x))
      waves <- numeric(length(id))
      for(i in 1:length(idIndex)) {
        waves[idIndex[[i]]] <- 1:length(idIndex[[i]])
      }
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
  dat$id <- id
  dat$weights <- weights
  dat$waves <- waves
  dat$N <- N
  dat$off <- off
  dat$sub.population <- sub.population
  subpopInd <- as.numeric(dat$sub.population)
  uniqueSubpop <- sort(unique(subpopInd))
  subpopInd <- sapply(subpopInd, function(x) which(x == uniqueSubpop))
  dat$subpopInd <- subpopInd

  # Checking that the treatment variables is in the model
  ind <- 0
  for(i in 1:(ncol(dat) + 1)) {
    if(i <= ncol(dat)) {
      if(all(treatment == dat[, i])) ind <- i
    } else if(ind == 0) {
      stop("treatment variable must be in the model!")
    }
  }
  treatmentColIndex <- ind

  dat <- dat[with(dat, order(id, waves)), ]

  # Creating 'full' dataset with all id/waves combinations
  imputedDat <- with(dat, expand.grid(id = unique(id), waves = unique(waves)))
  imputedDat <- merge(imputedDat, dat, all = TRUE)
  imputedDat <- imputedDat[with(imputedDat, order(id, waves)), ]
  treatmentColIndex <- which(names(imputedDat) == names(dat)[treatmentColIndex])

  #### Preparing new dataset with multiple instances per id for different effect levels
  # Checking that the treatment variable fulfills conditions
  if(!is.numeric(treatment) & !is.logical(treatment) & is.integer(treatment)) {
    stop("Treatment should be an integer or boolean vector where 0 means no treatment effect")
  }
  else if(any(round(treatment) != treatment)) {
    stop("Treatment should be an integer or boolean vector where 0 means no treatment effect")
  }

  if(all(treatment != 0)){
    warning("No 0's in treatment vector to serve as baseline effect! Convergence could be slow!")
  }

  #### Augmenting dataset
  dataByID <- with(imputedDat, split(imputedDat, id))
  dataByID <- lapply(dataByID, augment.dataset, treatment.levels, treatmentColIndex)
  augmentedData <- do.call("rbind", dataByID)
  augmentedData[, treatmentColIndex] <- factor(augmentedData[, treatmentColIndex])
  augmentedData$mixtureID <- factor(paste("mix", augmentedData$id, augmentedData$treatmentLevel, sep = ""))
  augmentedData <- augmentedData[with(augmentedData, order(id, mixtureID, waves)), ]

  #### Analysis preliminaries
  id <- augmentedData$mixtureID
  uniqueID <- unique(id)
  idIndex <- lapply(uniqueID,function(x) which(id == x))
  levelsPerObs <- augmentedData$treatmentLevel

  rho <- init.rho
  waves <- augmentedData$waves
  uniqueWaves <- unique(waves)
  nUniqueWaves <- length(unique(waves))
  for(i in 1:length(uniqueWaves)) {
    waves[which(waves == (sort(uniqueWaves)[i]))] <- i
  }
  augmentedData$waves <- waves

  # Structured Correlation
  structured <- diag(max(waves))

  X <- model.matrix(formula, augmentedData)

  mixtureWeights <- augmentedData$mixtureWeights
  rawWeights <- augmentedData$weights

  sub.population <- augmentedData$sub.population
  uniqueSubpop <- unique(sub.population)
  uniqueSubpop <- uniqueSubpop[!is.na(uniqueSubpop)]
  N <- augmentedData$N
  y <- model.response(model.frame(formula, augmentedData))
  notNAind <- !is.na(mixtureWeights)
  if(max(y) > 1) {
    prop <- y/N[notNAind]
  } else {
    prop <- y
  }
  # Replacing zero proportions with small proportions
  prop[prop == 0] <- min(prop[prop != 0] / 10^2, na.rm = TRUE)
  augmentedData$prop[notNAind] <- prop
  count <- prop * N[notNAind]
  off <- augmentedData$off

  converged <- FALSE
  subjectProbs <- matrix(nrow = length(unique(dat$id)), ncol = treatment.levels)
  levelProbs <- matrix(1, nrow = length(uniqueSubpop), ncol= treatment.levels + 1)
  deviance <- 10^19
  eta <- NULL
  subjectOffset <- off
  beta <- 0

  # Optimization Loop
  for(iter in 1:maxiter) {
    mixtureWeights[is.na(mixtureWeights)] <- 0
    tempID <- id[notNAind]
    glmerTime <- system.time(fit <- lme4::glmer(cbind(y, N[notNAind]) ~ X - 1
                 + (1 | tempID),
                 etastart = eta,
                 family = binomial, nAGQ = nAGQ,
                 weights = mixtureWeights[notNAind]))

    eta <- predict(fit)
    fixedEta <- predict(fit, re.form = NA)
    mu <- expit(eta)
    randomsd <- sqrt(summary(fit)$varcor[[1]][1])
    newbeta <- coef(fit)[[1]][1, -1]

    augmentedData$mu[notNAind] <- mu
    augmentedData$eta[notNAind] <- eta
    augmentedData$fixedEta[notNAind] <- fixedEta
    augmentedData$randomMean[notNAind] <- eta - fixedEta

    # Updating  level probabilities
    for(k in 1:length(uniqueSubpop)) {
      p <- sapply(1:treatment.levels,
                  function(x) sum(mixtureWeights[levelsPerObs==x & sub.population==uniqueSubpop[k]],
                                  na.rm=TRUE))
      p <- p / sum(p)
      levelProbs[k, 2:(treatment.levels + 1)] <- p
    }

    ## Updating Mixture Weights
    if(iter == 1){
      nsamp <- 15
    } else {
      nsamp <- ceiling(ifelse(weightTime[3] < glmerTime[3], nsamp * 1.25, nsamp * 0.75))
      nsamp <- max(nsamp, 15)
    }
    print(nsamp)
    weightTime <- system.time(newWeights <- by(augmentedData, augmentedData$id,
                         computeGmmWeights, randomsd, nsamp = nsamp, levelProbs))
    newWeights <- unlist(newWeights)
    if(iter == 1) {
      mixtureWeights <- newWeights
    } else {
      mixtureWeights <- mixtureWeights + (newWeights - mixtureWeights) / sqrt(iter - 1)
    }
    mixtureWeights <- unlist(mixtureWeights)
    mixtureWeights <- pmax(mixtureWeights, 10^-3)
    mixtureWeights <- mixtureWeights * rawWeights
    augmentedData$mixtureWeights <- mixtureWeights

    # Deciding whether to stop
    print(iter)
    print(levelProbs)
    print(round(newbeta, 3))
    if(iter < 4) {
      print(max(abs(beta - newbeta)))
    }

    beta <- newbeta
    if(iter == 1) {
      betaMat <- matrix(ncol = length(beta), nrow = maxiter)
    }
    betaMat[iter, ] <- as.numeric(beta)

    if(iter >= 4) {
      betaInd <- ceiling(iter/3):iter
      betaMean <- colMeans(betaMat[betaInd, ])
      betaSD <- apply(betaMat[betaInd, ], 2, sd) / sqrt(length(betaInd))
      diff <- max((betaSD / max(abs(betaMean))))
      print(round(diff, 3))
      if(max(diff) <= tol) {
        break
      }
    }
  }

  # Computing the subject probabilities
  mixtureWeights[mixtureWeights == 10^-3] <- 0
  subjectProbs <- by(mixtureWeights, augmentedData$id, function(x) unique(x))
  subjectProbs <- do.call("rbind", subjectProbs)

  # Preparing Level Probabilities or Output
  levelProbs <- data.frame(levelProbs)
  names(levelProbs) <- paste("level",c(0:treatment.levels),sep="_")
  rownames(levelProbs) <- unique(sub.population[notNAind])

  gmmfit <- fit
  fit <- list()
  fit$subject.probs <- subjectProbs
  fit$call <- match.call()
  fit$beta <- beta
  fit$iterations <- iter
  fit$converged <- converged

  mu <- by(expit(eta), augmentedData$id, function(x) (matrix(x, ncol = 2)))
  mu <- do.call("rbind", mu)
  mu <- data.frame(id = dat$id, population = dat$sub.population,
                   waves = dat$waves, nullMu = mu[, 1], altMu = mu[, 2])
  fit$mu <- mu
  fit$eta <- eta
  randomIntercept <- matrix(unique(round(eta - fixedEta, 8)), byrow = TRUE, ncol = 2)
  fit$randomSD <- randomsd
  fit$gmmfit <- gmmfit
  fit$randomEffectEst <- randomIntercept
  fit$X <- X
  return(fit)
}
