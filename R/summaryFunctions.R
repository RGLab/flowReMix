#' @name plot
#' @title plot a flowReMix fit object
#' @description Generate a plot of a flowReMix model fit. Various plot types are supported.
#' \itemize{
#' \item{subject}{The name of the subject variable, as an unquoted variable.}
#' \item{target}{unquoted variable name in the data corresponding to the outcome variable. Required.}
#' \item{varname}{character string used to label the outcome variable on the plot legend.}
#' \item{type}{character specifying the plot type. can be "scatter", "boxplot", "FDR", "ROC", "ising". Required.}
#' }
#' @param x The model fit of class \code{flowReMix} returned by the fitting function.
#' @param ... additional arguments, see description.
#' @usage
#'   \method{plot}{flowReMix}(x,...)
#' @export
plot.flowReMix <- function(x, target = NULL, varname = NULL,
                           subsets = NULL, ncol = 5,
                           palette = NULL, paletteRange = NULL,
                           type = c("ROC","scatter","boxplot","FDR","graph"), ...){
  mc = match.call()
  type=type[1] #need to do this otherwise you get warnings due to different lengths in an equality comparison when using the default passed vector.
  if(!is.null(mc$target)){
    target = mc$target
    target = enquo(target)
    subject_id = x$subject_id
    target = x$data %>% group_by(!!subject_id) %>% summarize(outcome=unique(!!target))%>%ungroup#%>%select(outcome)%>%unlist%>%factor
    target <- as.data.frame(target)
    post <- x$posteriors[, 1:2]
    if(is.factor(x$posteriors[, 1])) {
      target[, 1] <- factor(target[, 1], levels = levels(x$posteriors[, 1]))
    }
    target <- merge(post, target, sort = FALSE)
    target <- target[, 3]
    mc$target = target
  }

  if(type == "FDR") {
    table <- fdrTable(obj = x, target = target)
    mc[[1]] = as.name("plot")
    mc$obj = table
    mc$target = NULL
    return(eval(mc,envir = parent.frame()))
  } else if(type == "ROC") {
    # mc[[1]] = as.name("flowReMix:::plotROC")
    mc[[1]] = getFromNamespace("plotROC",ns = "flowReMix")
    mc$obj = mc$x
    mc$x =NULL
    return(eval(mc,envir = parent.frame()))
  } else if(type == "scatter") {
    figure <- plotScatter(x, subsets = subsets, target = target,
                          varname = varname, ncol = ncol,
                          colPalette = palette, paletteRange = paletteRange)
    return(figure)
    # mc[[1]] = as.name("flowReMix:::plotScatter")
    # mc[[1]] = getFromNamespace("plotScatter",ns = "flowReMix")
    # mc$obj = mc$x
    # mc$x =NULL
    return(eval(mc,envir = parent.frame()))
  } else if(type == "boxplot") {
    # mc[[1]] = as.name("flowReMix:::plotBoxplot")
    mc[[1]] = getFromNamespace("plotBoxplot",ns = "flowReMix")
    mc$obj = mc$x
    mc$x =NULL
    return(eval(mc,envir = parent.frame()))
  } else if(type == "graph") {
    # mc[[1]] = as.name("flowReMix:::plotRawGraph")
    mc[[1]] = getFromNamespace("plotRawGraph",ns = "flowReMix")
    mc$obj=mc$x
    mc$x=NULL
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
#' \itemize{
#' \item{target}{The name of the outcome variable as an unquoted variable.}
#' \item{type}{either "ROC" or "FDR".}
#' }
#' @param object a flowReMix model fit.
#' @param ... additional arguments. See description
#' @importFrom rlang enquo
#' @export
summary.flowReMix <- function(object, target, type = c("ROC", "FDR"), direction = "auto", adjust = "BH",
                              test = c("wilcoxon", "logistic", "ttest"),
                              sortAUC = FALSE, ...) {
  mc = match.call();
  if(is.null(mc$subject_id)){
    subject_id = object$subject_id
    subject_id = enquo(subject_id)
  }else{
    subject_id = mc$subject_id
    subject_id = enquo(subject_id)
  }
  if(is.null(mc$target)){
    stop("Please specify an argument for `target`. \n This should be the unquoted
         name of an outcome variable in the data. \n
         e.g.: summary(fit, target = outcome)")
  }

  type <- match.arg(type[1], c("FDR","ROC"))
  if(!exists("data",object)){
    stop("modify the fit object to contain the input data as element `fit$data`")
  }

  #Check of the target variable is valid
  target <- match.call()$target
  isvalid = object$data %>%
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

  target <- match.call()$target
  outcome = object$data %>%
    group_by(!!subject_id) %>%
    summarize(outcome=unique(!!target))
  #left join ensures order in posteriors is respected.
  outcome <- as.data.frame(outcome)
  object$posteriors[, 1] <- as.character(object$posteriors[, 1])
  outcome[, 1] <- as.character(outcome[, 1])
  outcome = suppressWarnings(left_join(object$posteriors,outcome, by = quo_name(subject_id)) %>% select(outcome) %>%unlist)
  if("ROC" == type) {
    table <- rocTable(object, outcome, direction = direction, adjust = adjust,
                         pvalue = test, sortAUC = sortAUC,...)
    return(table)
  } else if(type == "FDR") {
    return(fdrTable(object, ifelse(is.factor(outcome),outcome,factor(outcome))))
  } else {
    stop("Unknown summary method!")
  }
}

#' @title print.flowReMix
#' @name print.flowReMix
#' @description print a flowReMix object
#' @param x a flowReMix fit
#' @param ... additional arguments
#' @export
print.flowReMix = function(x, ...){
  cat("A flowReMix fit with \n")
  cat("\t",ncol(x$posteriors)-1," cell subsets\n")
  cat("\t",nrow(x$posteriors)," subjects\n")
}






