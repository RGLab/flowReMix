gibbsForParallel <- function(subjectData, isingCoefs,
                             covariance, nsamp, nSubsets, keepEach,
                             intSampSize, MHcoef, M, betaDispersion,
                             invcov) {
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
                                                as.integer(subjectData$preAssignment))

  unifVec <- runif(nsamp * nSubsets)
  eta <- subjectData$nullEta
  assignment <- assignmentMat[nrow(assignmentMat), ]
  responderSubset <- popInd %in% which(assignment == 1)
  randomEst <- as.numeric(estimatedRandomEffects[i, ])
  eta[responderSubset] <- subjectData$altEta[responderSubset]
  success <- 0
  attempt <- 0
  randomMat <- randomEffectCoordinateMH(y, N,
                                        as.integer(i), nsamp, nSubsets,
                                        MHcoef,
                                        as.vector(assignment),
                                        as.integer(popInd),
                                        as.numeric(eta),
                                        randomEst,
                                        as.numeric(condvar), covariance, invcov,
                                        unifVec,
                                        M, betaDispersion,
                                        keepEach)

  return(list(assign = assignmentMat, random = randomMat, s = success, a = attempt))
}


#' Fitting a Mixture of Mixed Effect Models for Binomial Data
#'
#' @description \code{flowReMix} fits a mixture of mixed effect models
#'   to binomial or over-dispersed binomial data. The package was specifically
#'   designed for analyzing flow-cytometry cell-count data but may be suitable
#'   for other purposes as well.
#'
#' @param formula an object of class \code{\link[stats]{formula}}. The response
#'   should be a matrix of two column matrix with first column containing the
#'   counts of the cell subsets of interest and the second column the difference
#'   between the reference count and the cell count.
#'
#' @param cell_type a factor vector identifying which cell type each row in the
#'   data set refers to.
#'
#' @param subject_id a vector identifying the subjects.
#'
#' @param data a data frame containing the variables in the model. It is
#'   advisable to include the \code{subject_id}, \code{cell_type} and
#'   \code{cluster_variable} variables in the data frame.
#'
#' @param cluster_variable a variable with respect to which clustering will be
#'   done. See description for more detail.
#'
#' @param cluster_assignment an optional matrix of known cluster assignments.
#'   Must include all subject/cell_type combinations. See description for more
#'   detail.
#'
#' @param weights an option vector of weights.
#'
#' @param covariance the method to be used for estimating the covariance
#'   structure of the random effects. \code{\link[PDSCE]{pdsoft.cv}} will be
#'   used by default.
#'
#' @param ising_model a method for estimating the Ising model.
#'   \code{\link[IsingFit]{IsingFit}} will be used by default.
#'
#' @param regression_method the regression method to be used. Default option is
#'   the \code{\link[stats]{glm}} function with family = "binomial".
#'
#' @param iterations the number of stochastic-EM itreations to perform.
#'
#' @param verbose whether to print information regrading the fitting process as
#'   the optimization algorithm runs.
#'
#' @param control an optional object of \code{\link{flowReMix_control}} class.
#'
#'
#' @details flowReMix fits a mixture of mixed effects regression models for
#'   binomial data. Accordingly, the response supplied in the \code{formula}
#'   must contain be a two column matrix the first column of which is the number
#'   of successes and the second column is the number of failiures. In the
#'   context of flow-cytomery count the left column would be the cell counts and
#'   the right columns the parent counts minus the cell count. The right side of
#'   the formula should include any number of fixed effects. For details on how
#'   the function processes the formula object see, for example, the
#'   documentation for the \code{\link[stats]{glm}} function.
#'
#'   The model fit by the function is a hierchical one, assuming the existence
#'   of subjects and one or more cell-types for each subject. the
#'   \code{subject_id} variable identifies different rows in the dataset as
#'   corresponding to measurements taken from specific subjects. The model
#'   assumes the existence of a random intercept for each \code{cell_type}.
#'
#'   The \code{cluster_variable} identifies which variable out of the covariates
#'   corresponds to the variable with respect to which clustering should be
#'   performed. The model assumes that the effect of cluster variable (and
#'   corresponding interactions) are either always zero or non-zero. For
#'   flow-cytomery experiments the \code{cluster_variable} will typically be an
#'   indicator for whether the stimulation introduced into the blood sample was
#'   an antigen or a control. A response status (zero or non-zero) is estimated
#'   for each subject/cell-type combination. The dependence between the
#'   cell-subsets is modeled via an Ising model.
#'
#'   \code{cluster_assignment} is an optional variable which allows the user to
#'   pre-specificy some known cluster assignments. For example, in vaccine
#'   studies we could expect all subjects who received a placebo treatment to be
#'   non-responders across all cell-subsets. This variable should be three
#'   column matrix, the first column of which should contain all unique values
#'   of \code{subject_id}, the second should column should contain all unique
#'   values of \code{cell_type} and in total the matrix should include all
#'   \code{subject_id} and \code{cell_type} combinations. The third column is an
#'   integer which takes the value 0 if the cell-type/subject combination is
#'   non-response, 1 if it is response and -1 if the response status is unknown
#'   and must be estimated.
#'
#'   The fitting algorithm uses one of three methods for estimating the
#'   covariance structure of the random effects. A diagonal covariance structure
#'   will be estimated if \code{covariance = "diagonal"}. A dense covariance
#'   structure will be estimated with no penalization will be estimated if
#'   \code{covariance = "dense"}. This may produce a singual covariance
#'   structure if the number of subjects is smaller than the number of
#'   cell-types. A sparse covariance matrix is estimated via the
#'   \code{\link[PDSCE]{pdsoft.cv}} function by default.
#'
#'   The ising model describing the dependence between the response/non-resposne
#'   status of the different cell-types can be estimated via three methods. If
#'   \code{ising_model} is set to \code{"none"} then an independence model is
#'   assumed. If the \code{ising_model} is set to \code{"dense"} then the ising
#'   model is estimated via a set of firth regressions
#'   (\code{\link[logistf]{logistf}}), one for each node in the graph. The
#'   default option is \code{"sparse"}, where the
#'   \code{\link[IsingFit]{IsingFit}} method is used.
#'
#'   \code{regression_method} specifies which function should be used for
#'   estimating the reqression coefficients conditionally on the values of the
#'   random effects and cluster assignments. If the default option
#'   \code{"binom"} is chosen then a binomial model is fit using the
#'   \code{\link[stats]{glm}} function. Otherwise, if \code{"betabinom"} option
#'   is selected then a beta-binomial regression model is estimated with the
#'   \code{\link[gamlss]{gamlss}} function. We recommend using the
#'   \code{"sparse"} method which uses the \code{\link[glmnet]{cv.glmnet}}
#'   procedure if the number of subjects is small and the number of predictors
#'   is large.
#'
#' @return \code{flowReMix} returns an object of class \code{flowReMix} which
#'   contains the following variables:
#'
#'   * \code{coefficients} a list, each component of which is a vector of
#'   regession coefficients corresponding to a single cell type.
#'
#'   * \code{posteriors} a matrix containing the posterior probabilities for
#'   response computed for each subject/cell-type combination.
#'
#'   * \code{levelProbs} a vector of the marginal estimated probabilities of
#'   response estiamted for each cell subset.
#'
#'   * \code{randomEffects} the estimated random effects for each
#'   subject/cell-type.
#'
#'   * \code{covariance} the estimated covariance structure for the random
#'   effects.
#'
#'   * \code{isingCov} the estimated `covariance` structure of the ising model.
#'
#'   * \code{isingfit} if \code{ising_model = "sparse"} then the object returned
#'   by the \code{\link[IsingFit]{IsingFit}} function. \code{NULL} otherwise.
#'
#'   * \code{dispersion} the over-dispersion estimated for each cell-subset. If
#'   regression method is not "betabinomial" then this will be a vector of large
#'   constants.
#'
#'   * \code{assignmentList} a list of matrices containing the posterior cluster
#'   assignemnt sampled for each subject at the last iteration of the stochastic
#'   EM algorithm.
#'
#'
#' @md
#' @export
flowReMixP <- function(formula,
                      subject_id,
                      cell_type = NULL,
                      cluster_variable,
                      data = parent.frame(),
                      cluster_assignment = NULL,
                      weights = NULL,
                      covariance = c("sparse", "dense", "diagonal"),
                      ising_model = c("sparse", "dense", "none"),
                      regression_method = c("binom", "betabinom", "sparse"),
                      iterations = 10, verbose = TRUE,
                      control = NULL,
                      parallel = FALSE, ncores = 2) {
  # Initializing parallel computation
  if(parallel) {
    snowfall::sfInit(parallel = TRUE, cpus = ncores)
  }

  # Getting control variables -------------------
  if(is.null(control)) {
    control <- flowReMix_control()
  } else if(class(control) != "flowReMix_control") {
    stop("`control' variable must be of `flowReMix_control' class!")
  }
  updateLag <- control$updateLag
  randomAssignProb <- control$randomAssignProb
  nsamp <- control$nsamp
  dataReplicates <- control$nPosteriors
  maxDispersion <- control$maxDispersion
  centerCovariance <- control$centerCovariance
  intSampSize <- control$intSampSize
  initMHcoef <- control$initMHcoef
  keepEach <- control$keepEach
  initMethod <- control$initMethod
  if(is.null(initMethod)) {
    if(regression_method == "sparse") {
      initMethod <- "sparse"
    } else {
      initMethod <- "binom"
    }
  }

  # Checking if inputs are valid --------------------------
  regressionMethod <- regression_method
  isingMethod <- ising_model
  covarianceMethod <- covariance
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
  if(nsamp < keepEach) {
    nsamp <- keepEach
    warning("Number of samples per MH iteration must be equal to or larger than `keepEach` variable!")
  }

  call <- as.list(match.call())
  if(is.null(call$cluster_variable)) {
    stop("Treatment variable must be specified!")
  }

  if(length(isingMethod) > 1) isingModel<- isingModel[1]
  if(!(isingMethod %in% c("sparse", "dense", "none"))) {
    stop("ising_model must be one of sparse, dense or none")
  }

  if(length(regressionMethod) > 1) regressionMethod <- regressionMethod[1]
  if(!(regressionMethod %in% c(c("binom", "betabinom", "sparse")))) stop("regression_method must be one of binom, betabinom or sparse")
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
  if(!(covarianceMethod %in% c("sparse", "dense", "diagonal"))) {
    stop("`covariance' must be one of sparse, dense or diagonal!")
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
  initialization <- lapply(dataByPopulation, initializeModel, formula = initFormula, method = initMethod)
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
    preAssignment$preAssignment <- rep(-1, nrow(preAssignment))
  } else {
    preAssignment <- cluster_assignment
    if(ncol(preAssignment) != 3) stop("preAssignment must have 3 columns: id, subset, assignment")
    if(any(!(preAssignment[, 3] %in% c(-1, 0, 1)))) stop("The third column of preAssignment must take values -1, 0 or 1.")
    preAssignment <- data.frame(preAssignment)
    names(preAssignment) <-  c("id", "subset", "preAssignment")
    if(dataReplicates > 1) {
      preAssignment <- do.call("rbind", lapply(1:dataReplicates, function(x) replicateDataset(preAssignment, x)))
    }
    if(nrow(preAssignment) != (nSubsets * nSubjects)) stop("preAssignment must have nSubjects X nSubsets rows.")
    if(any(!(preAssignment[, 1] %in% unique(dat$id)))) stop("The first column of Preassignment must correspond to the id variable.")
  }

  preAssignment <- preAssignment[order(preAssignment$id, preAssignment$subset), ]
  dat <- merge(dat, preAssignment,
               by.x = c("id", "sub.population"),
               by.y = c("id", "subset"), all.x = TRUE)

  # More preparations ---------------------------
  dat$tempTreatment <- dat$treatmentvar
  databyid <- by(dat, dat$id, function(x) x)
  ### ADD RANDOM EFFECT TO DAT HERE !!! ----------------------
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

  clusterDensities <- matrix(nrow = 2, ncol = intSampSize)
  isingCoefs <- matrix(0, ncol = nSubsets, nrow = nSubsets)
  if(betaDispersion) {
    M <- rep(10^4, nSubsets)
  } else {
    M <- rep(10^8, nSubsets)
  }

  flagEquation <- rep(0, nSubsets)

  # Starting analysis -------------------------------
  for(iter in 1:maxIter) {
    # Refitting Model with current random effects/assignments ---------------
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

    # Updating Prediction --------------------------
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

    # S-step --------------------------
    gibbsResult <- snowfall::sfLapply(databyid, gibbsForParallel, isingCoefs, covariance, nsamp,
                                      nSubsets, keepEach, intSampSize, MHcoef,
                                      M, betaDispersion, invcov)

    for(i in 1:nSubjects) {
      assignmetMat <- gibbsResult[[i]]
      iterPosteriors <- colMeans(assignmentMat)
      posteriors[i, ] <- posteriors[i, ] + (iterPosteriors - posteriors[i, ])/max(iter - updateLag, 1)
      clusterAssignments[i, ] <- assignmentMat[nrow(assignmentMat), ]
      assignmentList[[i]] <- assignmentMat
    }
    for(i in 1:nSubjects) {
      # MH sampler for random effects
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
    if(verbose) cat("\n")

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
      if(regression_method == "betabinom") print(c(M = M))
      print(round(rbind(MH = MHcoef, ratevec = ratevec), 3))
    }
  }

  # Processing posteriors -------------------
  posteriors <- data.frame(posteriors)
  uniqueIDs <- sapply(databyid, function(x) x$id[1])
  if(dataReplicates <= 1) {
    posteriors <- cbind(id = uniqueIDs, 1 - posteriors)
    names(posteriors) <- c(as.character(call$subject_id), names(dataByPopulation))
  } else {
    realIDs <- gsub("\\%%%.*", "", uniqueIDs)
    post <- by(posteriors, INDICES = realIDs, FUN = colMeans)
    postid <- names(post)
    posteriors <- data.frame(do.call("rbind", post))
    names(posteriors) <- names(dataByPopulation)
    posteriors <- cbind(id = postid, 1 - posteriors)
    names(posteriors)[1] <- as.character(call$subject_id)
  }

  # Processing random effects -----------
  estimatedRandomEffects <- data.frame(estimatedRandomEffects)
  names(estimatedRandomEffects) <- names(coefficientList)
  estimatedRandomEffects <- cbind(id = uniqueIDs, estimatedRandomEffects)

  # Preparing flowReMix object --------------------
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
  class(result) <- "flowReMix"
  if(parallel) sfStop()
  return(result)
  }

