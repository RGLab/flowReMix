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

#' @export
summary.flowReMix <- function(obj, ...) {
  args = list(...)
  target = args[["target"]]
  type = args[["type"]]
  type = match.arg(type,c("FDR","ROC"))
  # some more error checking of target
  type <- type[1]
  if("ROC" == type) {
    return(rocTable(obj, target, ...))
  } else if(type == "FDR") {
    return(fdrTable(obj, target))
  } else {
    stop("Unknown summary method!")
  }
}






