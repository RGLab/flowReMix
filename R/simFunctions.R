flowFitToSim <- function(fit) {
  if(!("flowReMix" %in% class(fit))) {
    stop("`fit' object must be of class `flowReMix'")
  }
  dat <- buildFlowFrame(obj$call, obj$data)
  formula <- editFlowFormulas(fit$call, FALSE)$initFormula
  X <- model.matrix(formula, dat)
  coef <- fit$coefficients
  dispersion <- fit$dispersion
  ising <- fit$isingCov
  covariance <- fit$covariance

  simobj <- list(X = X, coef = coef, dispersion = dispersion, ising = ising, cov = covariance)
  class(simobj) <- "flowReMix_sim"
  return(simobj)
}
