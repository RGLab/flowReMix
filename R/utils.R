#'@importFrom stats model.frame model.matrix na.omit na.pass quantile sd
#'@importFrom utils getFromNamespace
NULL

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

#' getIsing
#'
#' @param x \code{flowReMix} model fit.
#'
#' @return \code{isingStability} slot from the model.
#' @export
#'
#' @examples
#' data(fit505)
#' getIsing(fit505)
getIsing = function(x){
  .isFlowRemix(x)
  return(x$isingStability)
}


#' degreeFromStringFun
#' A function to parse the degree of a subset from its name.
#' @details Assumes that a cell subset name follows a specific pattern where
#' the positive functions are simply listed in the name. For example, "8+/A+B+C+" would be  a degree 3
#' subset of CD8+ cells.
#' @param x \code{character} cell subset name from \code{getSubsets}
#'
#' @return \code{numeric} vector of degree of functionality.
#' @export
#' @seealso \code{\link{flowReMixPFS}} \code{\link{weightForPFS}}
#' @examples
#' data(fit505)
#' degreeFromStringFun(getSubsets(fit505))
#'
degreeFromStringFun = function(x){
  degree = unlist(
    Map(
      strsplit(x = gsub(".*?/", "", x), split = "\\+"),
      f=length),
    use.names = FALSE)
  names(degree) = x
  degree
}

#' weightForPFS
#' Compute the weighting for the polyfunctionality score for each type of subset.
#' @details Computes the weight for the polyfunctionality score calculation. The weight is \code{choose(M,k)} where
#' k is the degree of a subset and M is the \code{max(degree)} or the number of markers.
#' @param x \code{flowReMix} fit object
#' @param M \code{numeric} the max possible functionality. Uses \code{max(degree)} if not provided.
#' @param parser \code{function} to parse the degree from a cell subset name. Defaults to \code{\link{degreeFromStringFun}} in the package.
#'
#' @return \code{numeric} vector of weights.
#' @seealso \code{\link{flowReMixPFS}} \code{\link{degreeFromStringFun}}
#' @export
#'
#' @examples
#' data(fit505)
#' weightForPFS(fit505)
weightForPFS = function(x, M = NULL, parser = degreeFromStringFun){
  if("character"%in%class(x)){
    degree =  parser(x)
  }else{
    .isFlowRemix(x)
    degree = parser(getSubsets(x))
  }
  if(is.null(match.call()$M))
     M = max(degree)
  return(degree/choose(M,degree))
}


#' perSubsetPFS
#' Compute the PFS (polyfunctionality score) per subject and stimulation group
#' @param x \code{flowReMix} model object
#' @param M \code{numeric} the max possible cell subset degree, the number of functions measured.
#' @param stimVar \code{name} Unquoted  name of the stimulation variable in the data. e.g. stim
#' @param parentVar \code{name} Unquoted name of the parent cell population variable in the data. e.g. parent
#' @param outcomeVar \code{name} Unquoted name of the outcome variable in the data. If provided, the scores will be merged with the data using the \code{subject_id}
#' @param ... additional arguments passed to weightForPFS.
#' @details
#' Requires that the data table has a variable for stimulation and for cell population parent.
#'
#' The user can pass in a function to the argument  \code{parser=}. This function parses a cell subset name string  and
#' returns the degree of functionality for the cell subset. The default is \code{\link{degreeFromStringFun}}
#' @note Uses rlang quosures to pass variable information. PFS are NOT normalized to the total number of possible subsets (n*(n+1)/2), like in COMPASS.
#' @seealso \code{\link{degreeFromStringFun}} \code{\link{weightForPFS}}
#' @return \code{data.frame} of polyfunctionality scores for each cell subset and subject, weighted appropriately.
#' @export
#'
#' @examples
#' data(fit505)
#' flowReMixPFS(fit505,M=5, stimVar = stimGroup, parentVar = parent)
#' flowReMixPFS(fit505,M=5,stimVar=stimGroup,parentVar=parent,outcomeVar=infection)
flowReMixPFS = function(x,M,stimVar = NULL, parentVar = NULL, outcomeVar = NULL, ...){
  mc  =  match.call()
  if(is.null(mc$stimVar)|is.null(mc$parentVar)){
    stop("stimVar and parentVar must be provided")
  }
  stimVar = enquo(stimVar)
  parentVar = enquo(parentVar)
  post = getPosteriors(x) %>% gather(pop,posterior,-1) %>% mutate(score=posterior*weightForPFS(pop,M=M, ...))
  subjid = x$subject_id
  subjid = enquo(subjid)
    post = x$data %>%
      select(!!parentVar,pop=subset,!!stimVar) %>%
      mutate(pop=as.character(pop)) %>%
      unique() %>%
      inner_join(post) %>%
      group_by(!!parentVar,!!stimVar, !!subjid) %>%
      summarize(PFS=mean(score))
    if(!is.null(mc$outcomeVar)){
      outcomeVar = enquo(outcomeVar)
      outcome = x$data%>%select(!!subjid,!!outcomeVar)%>%unique()
      post = inner_join(outcome,post)
    }

  return(post)
}


