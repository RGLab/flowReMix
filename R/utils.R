.isFlowRemix = function(x){
  if(!inherits(x,"flowReMix")){
    stop("x must be a flowReMix object.",call. = FALSE)
  }
}

#' getPosteriors
#'
#' @param x \code{flowReMix} model fit
#'
#' @return \code{data.frame} of posterior probabilities for each subject and subset in the model.
#' @export
#'
#' @examples
#' data(fit505)
#' getPosteriors(fit505)
#'
getPosteriors = function(x){
  .isFlowRemix(x)
  return(x$posteriors)
}

#' getSubsets
#'
#' @param x \code{flowReMix} object
#'
#' @return \code{vector} of subset names fitted in the model.
#' @export
#'
#' @examples
#' data(fit505)
#' getSubsets(fit505)
#'
getSubsets = function(x){
  .isFlowRemix(x)
  return(colnames(x$posteriors)[-1L])
}


