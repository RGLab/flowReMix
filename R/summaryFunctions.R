#' @name plot
#' @title plot a flowReMix fit object
#' @description Generate a plot of a flowReMix model fit. Various plot types are supported.
#' "scatter", "boxplot", "FDR", "ROC", "ising"
#' @param obj The model fit of class \code{flowReMix} returned by the fitting function.
#' @param type the type of plot to make, one of c("scatter","boxplot","FDR","ROC","graph")
#' @param target the name of the outcome variable as an unquoted variable.
#' @param varname the variable name to appear in the legend
#' @param ... additional arguments.
#' @export
plot.flowReMix <- function(obj,...){
  mc = match.call()
  if(!is.null(mc$target)){
    target = mc$target
    target = enquo(target)
    subject_id = obj$subject_id
    target = obj$data %>% group_by(!!subject_id) %>% summarize(outcome=unique(!!target))%>%ungroup%>%select(outcome)%>%unlist%>%factor
    mc$target = target
  }
  type = mc$type
  mc$type = NULL

  if(type == "FDR") {
    table <- fdrTable(obj, target = target)
    mc[[1]] = as.name("plot")
    mc$obj = table
    mc$target = NULL
    return(eval(mc,envir = parent.frame()))
  } else if(type == "ROC") {
    mc[[1]] = as.name("plotROC")
    return(eval(mc,envir = parent.frame()))
  } else if(type == "scatter") {
    mc[[1]] = as.name("plotScatter")
    return(eval(mc,envir = parent.frame()))
  } else if(type == "boxplot") {
    mc[[1]] = as.name("plotBoxplot")
    return(eval(mc,envir = parent.frame()))
  } else if(type == "graph") {
    mc[[1]] = as.name("plotRawGraph")
    mc$target=NULL
    if(!is.null(match.call()$fill) & is.null(match.call()$fillName)) {
      return(eval(mc,envir = parent.frame()))
    } else {
      return(eval(mc,envir = parent.frame()))
    }
  }
}


#' @name summary
#' @title summary of a flowReMix fit
#' @description summarize the output of a flowReMix object into a rocTable
#' Uses non-standard evaluation
#' @param subject the name of the subject variable, as an unquoted variable
#' @param target the name of the outcome variable as an unquoted variable.
#' @param type either "ROC" or "FDR".
#' @importFrom rlang enquo
#' @export
summary.flowReMix <- function(obj, ...) {
  mc = match.call();
  if(is.null(mc$subject_id)){
    subject_id = fit$subject_id
    subject_id = enquo(subject_id)
  }else{
    subject_id = mc$subject_id
    subject_id = enquo(subject_id)
  }
  if(is.null(mc$target)){
    stop("Please specify an argument for `target`. \n This should be the unquoted
         name of an outcome variable in the data. \n
         e.g.: summary(fit, target = outcome)")
  }else{
    target = mc$target
    target = enquo(target)
  }
  if(is.null(mc$type)){
    type = "ROC"
  }else{
    type = eval(mc$type,envir=parent.frame())
  }
  type = match.arg(type, c("FDR","ROC"))
  if(!exists("data",fit)){
    stop("modify the fit object to contain the input data as element `fit$data`")
  }
  #Check of the target variable is valid
  isvalid = obj$data %>%
    group_by(!!subject_id) %>% mutate(nlevels = length(unique(!!target)))
  if(!all(isvalid$nlevels %in% 1)){
    stop(
      quo_name(target),
      " must have one unique value per ",
      quo_name(subject_id),
      ". Found ",
      ifelse(
        length(unique(isvalid$nlevels)) == 1,
        unique(isvalid$nlevels),
        paste(range(unique(isvalid$nlevels)), collapse = " - ")
      )
    )
  }
  outcome = obj$data %>%
    group_by(!!subject_id) %>%
    summarize(outcome=unique(!!target))
  #left join ensures order in posteriors is respected.
  outcome <- as.data.frame(outcome)
  obj$posteriors[, 1] <- as.character(obj$posteriors[, 1])
  outcome[, 1] <- as.character(outcome[, 1])
  outcome = suppressWarnings(left_join(obj$posteriors,outcome, by = quo_name(subject_id)) %>% select(outcome) %>%unlist)
  if("ROC" == type) {
    mc$type = NULL
    mc[[1]]=as.name("rocTable")
    mc$target = outcome
    eval(mc, envir=parent.frame())
    # return(rocTable(obj, outcome, type=type,...))
  } else if(type == "FDR") {
    return(fdrTable(obj, ifelse(is.factor(outcome),outcome,factor(outcome))))
  } else {
    stop("Unknown summary method!")
  }
}







