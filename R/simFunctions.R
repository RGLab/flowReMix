#' @export
flowFitToSim <- function(fit, sampleSize = NULL, seed = NULL) {
  if(!("flowReMix" %in% class(fit))) {
    stop("`fit' object must be of class `flowReMix'!")
  }
  dat <- buildFlowFrame(fit$call, fit$data)$frame
  formula <- editFlowFormulas(fit$call, FALSE)$initFormula
  X <- model.matrix(formula, dat)
  coef <- fit$coefficients
  dispersion <- fit$dispersion
  ising <- fit$isingCov
  covariance <- fit$covariance
  call <- fit$call

  simobj <- list(X = X, coef = coef, dispersion = dispersion,
                 ising = ising, cov = covariance, formula = formula,
                 data = dat,
                 call = fit$call,
                 control = fit$control)
  class(simobj) <- "flowReMix_sim"
  if(!is.null(sampleSize)) {
    simobj <- generateFlowDataset(simobj, sampleSize, seed)
  }
  return(simobj)
}

#' @import IsingSampler
#' @import mvtnorm
#' @importFrom IsingSampler IsingSampler
#' @importFrom mvtnorm rmvnorm
#' @export
generateFlowDataset <- function(simobj, sampleSize, seed = NULL) {
  if(!is.null(seed)) {
    set.seed(seed)
  }

  if(!("flowReMix_sim" %in% class(simobj))) {
    stop("`fit' object must be of class `flowReMix_sim'!")
  }

  thresholds <- diag(simobj$ising)
  ising <- simobj$ising
  diag(ising) <- 0
  z <- IsingSampler(n = sampleSize, nIter = sampleSize * 20, graph = ising,
                    thresholds = thresholds, method = "MH")
  v <- rmvnorm(sampleSize, sigma = simobj$cov)
  data <- simobj$data
  data$id <- as.numeric(data$id)
  uniqueIds <- data$id
  dispersion <- simobj$dispersion

  if(is.factor(data$treatmentvar)) {
    baseline <- levels(data$treatmentvar)
  } else {
    baseline <- 0
  }

  subjectList <- list()
  for(i in 1:sampleSize) {
    sampId <- sample(uniqueIds, 1)
    subjDat <- subset(data, id == sampId)
    coef <- simobj$coefs[[i]]
    uniqueSubpop <- unique(subjDat$sub.population)
    poplist <- list()
    for(j in 1:length(uniqueSubpop)) {
      subset <- names(simobj$coef)[j]
      popDat <- subset(subjDat, sub.population == subset)
      if(z[i, j] == 1) {
        popDat$treatmentvar <- baseline
      }
      X <- model.matrix(simobj$formula, data = popDat)
      eta <- X %*% simobj$coef[[j]] + v[i, j]
      popDat$treatmentvar <- popDat$tempTreatment
      popDat$mu <- expit(eta)
      popDat$prob <- with(popDat, rbeta(length(mu), mu * dispersion[j], (1 - mu) * dispersion[j]))
      popDat$y <- with(popDat, rbinom(length(prob), N, prob))
      popDat$z <- z[i, j]
      popDat$v <- v[i, j]
      poplist[[j]] <- popDat
    }
    subjDat <- do.call("rbind", poplist)
    subjDat$id <- runif(1)
    subjectList[[i]] <- subjDat
  }
  simdata <- do.call("rbind", subjectList)
  simdata$prop <- simdata$y / simdata$N
  simdata$id <- factor(simdata$id)
  simobj$simdata <- simdata
  return(simobj)
}

#' @export
fitFlowSim <- function(simobj) {
  if(is.null(simobj$simdata)) {
    stop("Use `generateFlowDataset' to generate a sim object with a dataset!")
  }
  data <- simobj$simdata
  call <- simobj$call
  form <- update.formula(call$formula, cbind(y, N - y) ~ .)
  clustervar <- call$cluster_variable
  newfit <- flowReMix(formula = form,
                      subject_id = id,
                      cell_type = sub.population,
                      cluster_variable = clustervar,
                      data = data,
                      cluster_assignment = (call$cluster_assignment != FALSE),
                      covariance = call$covariance,
                      ising_model = call$ising_model,
                      regression_method = call$regression_method,
                      iterations = call$iterations,
                      parallel = call$parallel,
                      verbose = call$verbose,
                      control = simobj$control,
                      keepSamples = call$keepSamples)
  return(newfit)
}

