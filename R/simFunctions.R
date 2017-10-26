#' @export
flowFitToSim <- function(fit) {
  if(!("flowReMix" %in% class(fit))) {
    stop("`fit' object must be of class `flowReMix'!")
  }
  dat <- buildFlowFrame(obj$call, obj$data)
  formula <- editFlowFormulas(fit$call, FALSE)$initFormula
  X <- model.matrix(formula, dat)
  coef <- fit$coefficients
  dispersion <- fit$dispersion
  ising <- fit$isingCov
  covariance <- fit$covariance
  call <- fit$call

  simobj <- list(X = X, coef = coef, dispersion = dispersion,
                 ising = ising, cov = covariance, formula = formula,
                 data = dat)
  class(simobj) <- "flowReMix_sim"
  return(simobj)
}

#' @import IsingSampler
#' @import mvtnorm
#' @importFrom IsingSampler IsingSampler
#' @importFrom mvtnorm rmvnorm
#' @export
generateFlowDataset <- function(simobj, sampleSize) {
  if(!("flowReMix_sim" %in% class(simobj))) {
    stop("`fit' object must be of class `flowReMix_sim'!")
  }

  thresholds <- diag(simobj$ising)
  ising <- simobj$ising
  diag(ising) <- 0
  z <- IsingSampler(n = sampleSize, niter = sampleSize * 20, graph = ising,
                    thresholds = thresholds, method = "MH")
  v <- rmvnorm(sampSize, sigma = simobj$covariance)
  for(i in 1:length(simobj$coefs)) {
    coef <- simobj$coefs[[i]]
    subset <- names(simobj$coefs)[i]
    subX <- X[, ]
  }
}


