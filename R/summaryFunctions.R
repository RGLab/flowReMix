#' @export
plot.flowReMix <- function(obj, target = NULL, varname = NULL,
                           type = c("scatter", "boxplot", "FDR", "ROC", "ising"), ...) {
  if(!is.null(target) & is.null(varname)) {
    varname <- as.character(match.call()$target)
    varname <- varname[length(varname)]
  }

  if(is.null(target)) {
    target <- rep(1, nrow(obj$posteriors))
  }

  type <- type[1]
  if(type == "FDR") {
    table <- fdrTable(obj, target)
    return(plot(table, target = target, varname = varname, ...))
  } else if(type == "ROC") {
    return(plotROC(obj, target = target, varname = varname, ...))
  } else if(type == "scatter") {
    return(plotScatter(obj, target = target, varname = varname, ...))
  } else if(type == "boxplot") {
    return(plotBoxplot(obj, target = target, varname = varname, ...))
  } else if(type == "graph") {
    if(!is.null(match.call()$fill) & is.null(match.call()$fillName)) {
      fillName <- as.character(match.call()$fill)
      fillName <- fillName[length(fillName)]
      return(plotRawGraph(obj, fillName = fillName, ...))
    } else {
      return(plotRawGraph(obj, ...))
    }
  }
}


#' @name summary
#' @title summary of a flowReMix fit
#' @description summarize the output of a flowReMix object into a rocTable
#' Uses non-standard evaluation
#' @param subject the name of the subject variable, as an unquoted variable
#' @param target the name of the outcome variable as an unquoted variable. Default is outcome.
#'
#' @export
summary.flowReMix <- function(obj, subject = ptid, target = outcome,type="ROC",...) {
  type = match.arg(type,c("FDR","ROC"))
  target = enquo(target)
  subject = enquo(subject)
  if(!exists("data",fit)){
    stop("modify the fit object to contain the input data as element `fit$data`")
  }
  outcome = obj$data %>% group_by(!!subject) %>% summarize(outcome=unique(!!target))
  #ensure the order of outcome is as in posteriors.
  #NOTE Code should be changed to leverage merge etc throughout and not rely on rownames.
  outcome = as.data.frame(outcome)
  rownames(outcome) = outcome[,1]
  outcome = outcome[obj$posteriors[,1],]
  # some more error checking of target
  type <- type[1]
  if("ROC" == type) {
    return(rocTable(obj, outcome$outcome, ...))
  } else if(type == "FDR") {
    return(fdrTable(obj, outcome$outcome))
  } else {
    stop("Unknown summary method!")
  }
}






