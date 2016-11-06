wtQuantile <- function(quantile, x, wt) {
  wt <- wt / sum(wt)
  wt <- wt[order(x)]
  x <- sort(x)
  wt <- cumsum(wt)
  x <- c(x[max(which(wt <= quantile))], x[min(which(wt >= quantile))])
  return(mean(x))
}

generate.Ipp <- function(p) {
  I <- matrix(0,nrow=p^2,ncol=p^2)
  indices <- 1:p + p*(1:p-1)
  indices <- as.matrix(expand.grid(indices,indices))
  I[indices] <- 1
  return(I)
}

# A function for computing the t degrees of freedom
compute.t.df <- function(gradients, Ainv, weights, shrinkage.weight,
                         shrink.sandwich, shrunkSandwich = NULL) {
  vecB <- apply(gradients, 2, function(x) as.vector(x %*% t(x)))
  meanVecB <- rowSums(vecB)
  vecB <- apply(vecB, 2, function(x) x-meanVecB)
  covVecB <- 0
  for(i in 1:ncol(vecB)) covVecB <- covVecB + vecB[, i] %*% t(vecB[, i])
  covVecB <- covVecB / (sum(weights) * (sum(weights) - 1))
  covSandwich <- kronecker(Ainv, Ainv)
  covSandwich <- sum(weights)^2 * covSandwich %*% covVecB %*% covSandwich
  tdf <- diag(covSandwich)[seq(from=1,to=ncol(covSandwich),by=nrow(shrunkSandwich))+ (0:(nrow(shrunkSandwich)-1))]
  if(shrink.sandwich) {
    tdf <- 2*diag(shrunkSandwich)^2/(tdf*shrinkage.weight^2)
  } else {
    tdf <- 2*diag(shrunkSandwich)^2/tdf
  }

  #computing nu
  p <- nrow(shrunkSandwich)
  omegahat <- (diag(p^2) + generate.Ipp(p)) %*% kronecker(shrunkSandwich,shrunkSandwich)
  omegahat <- as.vector(omegahat)
  vecCovSandwich <- as.vector(covSandwich)#*shrinkage.weight^2
  nu <- t(vecCovSandwich)%*%omegahat/(t(vecCovSandwich)%*%vecCovSandwich)

  return(list(tdf = tdf, nu = nu))
}

# subject effect Stein shrinkage
stein.shrinkage <- function(x, weights, robust.param) {
  temp <- x
  signx <- sign(x)
  x <- abs(x)
  quantiles <- quantile(x, probs = c(1 - robust.param))
  x <- pmin(x, quantiles[1])
  x <- x * signx
  variance <- cov.wt(matrix(x, ncol = 1), weights)
  variance <- as.numeric(variance$cov)
  steinCoef <- max(1 - ((length(x) - 2) * variance / sum(x)^2), 0)
  return(x * steinCoef)

#   orderx <- order(x)
#   quantile <- sapply(c(robust.param, 0.5, 1 - robust.param), wtQuantile, x, weights)
#   x <- pmin(quantile[3], pmax(x, quantile[1]))
#   robustvar <- as.numeric(cov.wt(matrix(x, ncol = 1), weights)$cov)
#   zquant <- orderx / sum(!is.na(x))
#   x <- qnorm(zquant, sd = sqrt(robustvar))
  return(x)
}

# A function for computing the amount of shrinakge to perform on sandwich matrix
estimate.shrinkage.coefficient <- function(residMat, id,
                                           structured, weightMat,
                                           nfolds=5) {
  compute.cv.score <- function(lambda, whichFolds, foldResiduals,
                               unstructuredCorrelations, structured) {
    cvScores <- numeric(length(whichFolds))
    for(g in 1:length(whichFolds)) {
      Kg <- length(whichFolds)
      Rg <- as.matrix(lambda * unstructuredCorrelations[[g]] + (1 - lambda) * structured)
      logDet <- determinant(Rg,log=TRUE)[[1]][[1]]
      error <- sum(diag((foldResiduals[[g]] %*% solve(Rg) %*% t(foldResiduals[[g]]))), na.rm=TRUE)

      cvScores[g] <- -Kg/2*logDet -0.5*error
    }
    return(sum(cvScores))
  }

  noNA <- !is.na(weightMat)
  residMat <- residMat[noNA, ]
  weightMat <- weightMat[noNA]

  folds <- as.vector(sapply(1:nfolds, function(x)
    rep(x, ceiling((nrow(residMat) / nfolds)))))
  folds <- sample(folds, nrow(residMat), replace=FALSE)
  whichFolds <- lapply(1:nfolds,function(x) which(folds==x))

  unstructuredCorrelations <- lapply(whichFolds,
                                     function(x, residMat, weightMat) {
                                       cov.wt(residMat[-x, ], wt = weightMat[-x], cor = TRUE)$cor
                                     }, residMat, weightMat)

  foldResiduals <- lapply(whichFolds,function(x) residMat[x,])

  weight <- optimize(f=compute.cv.score,interval=c(0,.995),maximum=TRUE,
                     whichFolds=whichFolds,foldResiduals=foldResiduals,
                     unstructuredCorrelations=unstructuredCorrelations,
                     structured=structured)$maximum

  return(weight)
}

# A function for creating an augmented dataset for mixture modeling.
augment.dataset <- function(subjectData, treatment.levels, treatmentColIndex) {
  # if subject has no treatment do nothing
  notNA <- !is.na(subjectData[, ncol(subjectData)])
  if(all(subjectData[notNA, treatmentColIndex]==0)) {
    subjectData$treatmentLevel <- 0
    subjectData$mixtureWeights <- subjectData$weights
    return(subjectData)
  }

  subjectData <- lapply(1:treatment.levels, function(x) subjectData)
  initWeights <- runif(treatment.levels, min = 0.3, max = 0.7)
  initWeights <- initWeights / sum(initWeights)
  for(i in 1:treatment.levels) {
    if(i==1) subjectData[[i]][, treatmentColIndex] <- 0
    if(i>2) {
      nonZero <- subjectData[[i]][, treatmentColIndex] != 0
      subjectData[[i]][nonZero, treatmentColIndex] <- subjectData[[i]][nonZero, treatmentColIndex] + 0.01*i
    }
    subjectData[[i]]$mixtureWeights <- subjectData[[i]]$weights * initWeights[i]
    subjectData[[i]]$treatmentLevel <- i
  }

  subjectData <- do.call("rbind", subjectData)
  return(subjectData)
}

#' logit and expit functions
#'
#' @name logit
#' @aliases expit
#' @rdname logitexpit
#' @export
logit <- function(x) log(x / (1 - x))

#' @export
#' @rdname logitexpit
expit <- function(x) 1/(1 + exp(-x))

# a function for computing a global dispersion
compute.global.rho <- function(prop, mu, weights, N, robust.param, robust = TRUE) {
  presidsumsq <- function(rho) {
    var <- mu * (1 - mu) * (1 / N + (N - 1) / N * rho)
    presid <- residsq / var
    result <- weighted.mean(presid, weights, na.rm = TRUE) - 1
    #print(c(rho, result))
    return(result)
  }
  residsq <- (prop - mu)^2
  if(robust) {
    threshold <- sapply(c(robust.param/2, 1 - robust.param/2), wtQuantile, residsq, weights)
    residsq <- pmin(threshold[2], pmax(threshold[1], residsq))
  }
  result <- NULL
  try(result <- uniroot(f = presidsumsq, interval = c(10^-6, 0.5), tol = 10^-12)$root,
      silent = TRUE)
  if(is.null(result)) result <- 10^-6
  return(result)
}

# A function for computing subject level offsets
compute.subject.offset <- function(data, nWaves) {
  robust.param <- 0.001
  eta <- data$eta
  empeta <- with(data, logit(prop))
  eta <- matrix(eta, ncol = nWaves, byrow = TRUE)
  empeta <- matrix(empeta, ncol = nWaves, byrow = TRUE)
  randomEst <- empeta - eta
  randomEst <- rowMeans(empeta - eta, na.rm = TRUE)
  randomEst <- sapply(randomEst, function(x) rep(x, nWaves))
  return(as.vector(randomEst))
}

# A function for recalibrating pvalues.
calibrate.pvals <- function(nullPvals, pvals) {
  logPvals <- (pmax(pvals, 10^-5))
  logNullPvals <- (pmax(nullPvals, 10^-5))

  cdf <- ks::kcde(logNullPvals)
  whichMin <- sapply(logPvals, function(x) which.min(abs(x - cdf$eval.points)))
  recalibrated <- cdf$estimate[whichMin]
  #recalibrated <- sapply(logPvals, function(x) sum(logNullPvals <= x)) / length(nullPvals)

  return(recalibrated)
}

# copmute the mixture weighs from the likelihoods
compute.mixture.weights <- function(subjectData,
                                    nWaves,
                                    output = "weights") {
  likelihoods <- subjectData$likelihood
  likMat <- matrix(likelihoods, ncol = nWaves)
  likMat[is.na(likMat)] <- -Inf
  weights <- t(apply(likMat, 1, function(x) exp(x - min(x, na.rm = TRUE))))
  weights <- t(apply(weights, 1, function(x) x / sum(x, na.rm = TRUE)))
  weights[is.nan(weights)] <- NA

  if(output == "weights") {
    weights <- as.vector(weights)
  } else if(output == "subjectProbs") {
    weights <- (weights)
  } else {
    weights <- as.vector(apply(weights, 1, function(x) sum(log(x))))
    weights <- exp(weights - min(weights))
    weights <- weights / sum(weights)
  }
  weights[is.na(weights)] <- 1
  #weights <- pmax(weights, 0.1)
  return(weights)
}

# compute the non-full information gradient
compute.missing.information.grad <- function(index, weights, grads,
                                             levelsPerObs, sub.population,
                                             lambda, levelProbs) {
  subpop <- sub.population[index]
  uniqueSubpop <- unique(sub.population)
  uniqueSubpop <- uniqueSubpop[!is.na(uniqueSubpop)]

  weights <- weights[index]
  emgrads <- grads[, index]
  subjectLevels <- levelsPerObs[index]
  nlevels <- ncol(levelProbs) - 1
  nsubpop <- nrow(levelProbs)

  if(length(unique(subjectLevels)) == 1) {
    grads <- rbind(matrix(0, nrow=nlevels * nsubpop,
                          ncol = ncol(emgrads)), emgrads)
    return(grads)
  }

  obsPerLevel <- length(subjectLevels) / length(unique(subjectLevels))
  grads <- matrix(0,nrow=nrow(emgrads) + nlevels * nsubpop,
                  ncol=obsPerLevel)
  for(i in 1:obsPerLevel) {
    ind <- seq(from=i,
               to = length(subjectLevels) - obsPerLevel + i,
               by=obsPerLevel)

    probgrad <- rep(0,length(levelProbs[, -1]))
    whichSubpop <- which(subpop[ind][1] == uniqueSubpop)
    subind <- c(2*whichSubpop - 1, 2*whichSubpop)
    probgrad[subind] <- weights[ind] / levelProbs[whichSubpop, subjectLevels[ind] + 1] -
      lambda[whichSubpop]

    grads[, i] <- c(probgrad,rowSums(emgrads[, ind]))
  }

  return(grads)
}

# copmute the sandwich estimator for the mixture
compute.sandwich.mixturegee <- function(prop, mu, M, X, weights,
                                        originalID,
                                        levelsPerObs,
                                        sub.population, levelProbs,
                                        shrink.sandwich = TRUE,
                                        shrinkage.weight = 1, computeT = FALSE) {
  InvLinkDeriv <- function(mu) {
    eta <- log(mu/(1-mu))
    exp(-eta)/(1+exp(-eta))^2
  }

  resid <- prop-mu
  var <- mu*(1-mu)/(1+M)
  resid <- resid/var

  notNA <- !is.na(resid)
  inverseLinkDerivative <- InvLinkDeriv(mu[notNA])
  resid <- resid[notNA]
  weights <- weights[notNA]
  grads <- sapply(1:nrow(X), function(i) resid[i]*X[i, ]* weights[i] * inverseLinkDerivative[i])

  uniqueSubpop <- unique(sub.population[notNA])
  lambda <- sapply(uniqueSubpop,function(x) sum(weights[x==sub.population], na.rm = TRUE))
  originalIDindex <- (1:length(originalID[notNA]))
  originalID <- originalID[notNA]
  originalIDindex <- by(originalIDindex, originalID, FUN = function(x) x)
  missGrads <- lapply(originalIDindex, compute.missing.information.grad, weights,
                      grads, levelsPerObs, sub.population, lambda, levelProbs)

  A <- do.call("cbind",missGrads)
  A <- A[-seq(from=1,to=nrow(levelProbs)*2-1,by=2),]
  Amean <- rowMeans(A)
  #A <- apply(A,2,function(x) x-Amean)
  A <- A %*% t(A)

  B <- lapply(missGrads,function(x) rowSums(x))
  B <- do.call("cbind",B)
  B <- B[-seq(from=1,to=nrow(levelProbs)*2-1,by=2),]
  Bmat <- B %*% t(B)

  Ainv <- solve(A)
  sandwich <- Ainv %*% Bmat %*% Ainv

  if(computeT) {
    shrunkSandwich <- shrinkage.weight * sandwich + (1 - shrinkage.weight) * Ainv
    tdf <- compute.t.df(B, Ainv, weights, shrinkage.weight,
                        shrink.sandwich, shrunkSandwich)
    nu <- tdf$nu
    tdf <- tdf$tdf
    tdf <- pmax(tdf,1)
  } else {
    tdf <- 10^4
    nu <- 10^4
  }

  return(list(A = A, B = B, sandwich = sandwich, tdf = tdf, nu = nu))
}

# Estimate a composite mixture GEE
compositeMixture <- function(formula, sub.population = NULL,
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
                             robust.param = 0.2,
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
    # M Step - regression coefficients
    if(all(rho == 0)) {
      regWeights <- mixtureWeights
    } else {
      regWeights <- mixtureWeights / rho
    }

    mixtureWeights[is.na(mixtureWeights)] <- 0
    if(iter == 1 & !is.null(init.beta)) {
      newbeta <- 10^6
      beta <- init.beta
      if(length(beta) != ncol(X)) {
        warning("Illegal init.beta value!")
        beta <- NULL
        next
      }
    } else {
#       fit <- try(glm(prop ~ X - 1 + offset(off[notNAind]),
#                      weights = regWeights[notNAind],
#                      etastart = eta[notNAind],
#                      family = "binomial",
#                      control = glm.control(maxit = 5, epsilon = 10^-4)),
#                  silent = TRUE)

      fit <- glm.fit(x = X, y = prop[notNAind], weights = mixtureWeights[notNAind],
              etastart = eta[notNAind],
              family = binomial(),
              intercept = FALSE,
              control = glm.control(maxit = 3, epsilon = 10^-4))

      newbeta <- fit$coefficients
      if(is.null(beta)) {
        beta <- newbeta
      } else {
        beta <- beta + (newbeta - beta) / sqrt(iter)
      }
    }

    eta <- X %*% beta
    mu <- expit(eta)
    augmentedData$mu[notNAind] <- mu
    augmentedData$eta[notNAind] <- eta

    # Estimating subject offset
    if(compute.offset) {
      subjectOffset <- compute.subject.offset(augmentedData, nWaves)
      subjectOffset <- stein.shrinkage(subjectOffset, mixtureWeights, robust.param)
      subjectOffset <- subjectOffset + off
    } else {
      subjectOffset <- off
    }
    augmentedData$subjectOffset <- subjectOffset
    mu <- expit(eta + subjectOffset[notNAind])
    augmentedData$mu[notNAind] <- mu

    # Estimating Dispersion
    resid <- (prop-mu)
    rho <- compute.global.rho(prop, mu, mixtureWeights, N, robust.param)
    newM <- (1-rho)/rho
    M <- M + (newM - M) / sqrt(iter)

    # Making sure dispersion is not too small
    augmentedData$M[notNAind] <- M

    # Updating  level probabilities
    for(k in 1:length(uniqueSubpop)) {
      p <- sapply(1:treatment.levels,
                  function(x) sum(mixtureWeights[levelsPerObs==x & sub.population==uniqueSubpop[k]],
                                  na.rm=TRUE))
      p <- p / sum(p)
      levelProbs[k, 2:(treatment.levels + 1)] <- p
    }

    ## Updating Mixture Weights
    likelihoods <- with(augmentedData,
                        dbeta(prop, shape1 = M * mu, shape2 = M * (1-mu), log = TRUE) +
                          dbinom(count, N, prop, log = TRUE))

    # Checking for infinite and zero likelihoods
    likelihoods[likelihoods == Inf] <- log(max(exp(likelihoods[likelihoods != Inf])) * 10)
    likelihoods[likelihoods == -Inf] <- log(min(exp(likelihoods[likelihoods > 0])) / 10)

    # Multiplying the likelihoods by the overall level probabilities
    for(k in 1:length(uniqueSubpop)) {
      ind <- which(sub.population == uniqueSubpop[k])
      likelihoods[ind] <- likelihoods[ind] + log(levelProbs[k, levelsPerObs[ind] + 1])
    }

    augmentedData$likelihood <- likelihoods

    # Mixture Weights
    mixtureWeights <- by(augmentedData, augmentedData$id, compute.mixture.weights, nWaves = nWaves, output = "weights")
    mixtureWeights <- unlist(mixtureWeights)
    mixtureWeights <- mixtureWeights * rawWeights
    mixtureWeights <- pmax(mixtureWeights, 10^-6)
    mixtureWeights <- pmin(mixtureWeights, 1)
    augmentedData$mixtureWeights <- mixtureWeights

    # Deciding whether to stop
    print(iter)
    print(M)
    print(beta)
    print(levelProbs)
    print(max(abs(beta - newbeta)))
    if(max(abs(beta - newbeta))< tol & iter > 2) {
      break
    }
  }

  # Computing the subject probabilities
  subjectProbs <- by(augmentedData, augmentedData$id,
                     compute.mixture.weights, nWaves = nWaves,
                     output = "subjectProbs")
  for(i in length(subjectProbs):1) {
    if(all(subjectProbs[[i]]==1, na.rm = TRUE))  {
      subjectProbs[[i]] <- matrix(0,  nrow=length(subjectProbs[[i]]), ncol=treatment.levels)
      subjectProbs[[i]][, 1] <- 1
    }
  }
  subjectProbs <- do.call("rbind", subjectProbs)
  subjectProbs <- cbind(imputedDat$id, subjectProbs)

  # Recomputing dispersion without subject level intercept
  rho <- with(augmentedData, compute.global.rho(prop, expit(eta), mixtureWeights, N, robust.param, robust = FALSE))
  #print(M)
  M <- (1-rho)/rho
  #print(M)
  augmentedData$M <- M

  # Computing weight for shrinking the Sandwich
  if(shrink.sandwich) {
    rho <- 1/(1 + augmentedData$M)
    mu <- expit(augmentedData$eta) #augmentedData$mu
    prop <- augmentedData$prop
    var <- (mu * (1-mu)) * rho
    presid <- (prop-mu)^2/var
    residMat <- matrix(presid, ncol = nWaves, byrow = TRUE)
    shrinkWeights <- unlist(by(augmentedData, augmentedData$id,
                               compute.mixture.weights,
                               nWaves = nWaves,
                               output = "shrinkageWeights"))
    shrinkage.weight <- mean(replicate(5,
                             estimate.shrinkage.coefficient(residMat, id, structured,
                                                            shrinkWeights, nfolds=shrinkage.folds)))
  } else {
    shrinkage.weight <- 0.9
  }

  # Computing the Sandwich estimate
  sandwich <- with(augmentedData,
                   compute.sandwich.mixturegee(prop, expit(eta), M,
                                               X, mixtureWeights,
                                               id,
                                               treatmentLevel,
                                               sub.population,
                                               levelProbs,
                                               shrink.sandwich,
                                               shrinkage.weight,
                                               computeT = computeT))
  A <- sandwich$A
  tdf <- sandwich$tdf
  if(computeT)   tdf <- tdf[-(1:length(levelProbs[,1]))]
  nu <- sandwich$nu
  sandwich <- sandwich$sandwich

  # Computing coefficient table
  if(shrink.sandwich | TRUE) {
    robust.var <- sandwich*shrinkage.weight + solve(A) * (1-shrinkage.weight)
  } else {
    robust.var <- sandwich
  }

  coefTable <- matrix(nrow=length(beta),ncol=4)
  coefTable[, 1] <- beta
  coefTable[, 2] <- sqrt(diag(robust.var)[-(1:length(levelProbs[,1]))])
  coefTable[, 3] <- coefTable[,1]/coefTable[,2]
  coefTable[, 4] <- 2*pt(-abs(coefTable[,3]), tdf)
  coefTable <- data.frame(coefTable)
  rownames(coefTable) <- names(beta)
  names(coefTable) <- c("Coefficient","SD","t-stat","p-val")

  # Preparing Level Probabilities or Output
  levelProbs <- data.frame(levelProbs)
  names(levelProbs) <- paste("level",c(0:treatment.levels),sep="_")
  rownames(levelProbs) <- unique(sub.population[notNAind])

  # Computing p-values
  pvals <- NULL
  computePvals <- TRUE
  if(computePvals) {
    tempMu <- mu[levelsPerObs <= 1]
    tempDisp <- augmentedData$M[levelsPerObs <= 1]
    tempProp <- prop[levelsPerObs <= 1]

    pvals <- 1-pbeta(tempProp, tempDisp * tempMu, tempDisp* (1-tempMu))
  }

  fit <- list()
  fit$coefTable <- coefTable
  fit$subject.probs <- subjectProbs
  fit$call <- match.call()
  fit$beta <- beta
  fit$iterations <- iter
  fit$converged <- converged
  fit$naiv.var <- solve(A)
  fit$sandwich <- sandwich
  fit$var <- robust.var
  fit$rho <- 1/(M+1)
  fit$mu <- mu
  fit$eta <- eta
  fit$shrinkage.weight <- shrinkage.weight
  fit$tdf <- tdf
  fit$nu <- fit$nu
  fit$pvals <- pvals
  subjectOffset <- subjectOffset - weighted.mean(subjectOffset, mixtureWeights)
  fit$randomSD <- sqrt(weighted.mean((subjectOffset)^2, mixtureWeights))
  #fit$deviance <- sum(deviance)
  fit$level.probabilities <- levelProbs

  subjectOffsetMat <- do.call("rbind",
                              by(subjectOffset, augmentedData$id, function(x) matrix(x, ncol = 2)))
  propMat <- do.call("rbind",
                              by(logit(prop), augmentedData$id, function(x) matrix(x, ncol = 2)))
  etaMat <- do.call("rbind",
                    by(augmentedData$eta, augmentedData$id, function(x) matrix(x, ncol = 2)))
  diagnostic <- round(cbind(subjectProbs, subjectOffsetMat, propMat, etaMat), 3)
  diagnostic <- data.frame(diagnostic)
  names(diagnostic) <- c("ptid", "prob0", "prob1", "offset0", "offset1", "empeta0", "empeta1", "eta0", "eta1")
  return(fit)
}
