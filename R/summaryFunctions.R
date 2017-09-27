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
#' @param type either "ROC" or "FDR".
#' @export
summary.flowReMix <- function(obj, ...) {
  arglist = quos(...)
  if(!("subject_id"%in%names(arglist))){
    subject_id = fit$subject_id
  }else{
    subject_id = arglist[["subject_id"]]
  }
  if(!("target"%in%names(arglist))){
    target = quo(outcome)
  }else{
    target = arglist[["target"]]
  }
  if(!("type"%in%names(arglist))){
    type = "ROC"
  }else{
    type = quo_name(arglist[["type"]])
  }
  type = match.arg(type,c("FDR","ROC"))
  # target = enquo(target)
  # subject_id = enquo(subject_id)
  if(!exists("data",fit)){
    stop("modify the fit object to contain the input data as element `fit$data`")
  }
  outcome = obj$data %>%
    group_by(!!subject_id) %>%
    summarize(outcome=unique(!!target))
  #left join ensures order in posteriors is respected.
  outcome = suppressWarnings(left_join(obj$posteriors,outcome, by = quo_name(subject_id)) %>% select(outcome) %>%unlist)
  if("ROC" == type) {
    return(rocTable(obj, outcome, type=type))
  } else if(type == "FDR") {
    return(fdrTable(obj, ifelse(is.factor(outcome),outcome,factor(outcome))))
  } else {
    stop("Unknown summary method!")
  }
}







