#' @import ggplot2
#' @import dplyr
#' @importFrom pROC roc
plotROC <- function(obj, target, direction = "auto",
                    ncols = 5,
                    thresholdPalette = NULL,
                    paletteRange = NULL,
                    subsets = NULL, varname = NULL) {
  notNA <- !is.na(target)
  post <- obj$posteriors[notNA, -1]
  if(!is.null(subsets)) {
    post <- post[, colnames(post) %in% subsets]
  }
  target <- target[notNA]
  uniqVals <- length(unique(target))
  if(length(unique(target)) != 2) {
    stop("Target variable must have exactly two values!")
  }

  if(is.null(thresholdPalette)) {
    thresholdPalette <- viridis::plasma(256)
  }

  p <- ncol(post)
  plotList <- list()
  aucs <- numeric(p)
  for(i in 1:p) {
    rocfit <- roc(target ~ post[, i], smooth = FALSE,
                        direction = direction)
    aucs[i] <- rocfit$auc
    fpr <-  1 - rocfit$specificities
    fprvals <- unique(fpr)
    temp <- data.frame(sens = rocfit$sensitivities, fpr = 1 - rocfit$specificities,
                          threshold = rocfit$thresholds)

    temp$threshold[temp$threshold < 0] <- 0
    temp$threshold[temp$threshold > 1] <- 1
    temp$subset <- paste(colnames(post)[i], ": ", round(aucs[i], 3), sep = "")
    plotList[[i]] <- temp
  }

  if(is.null(varname)) {
    varname <- as.character(match.call()$target)
    varname <- varname[length(varname)]
  }

  plotdat <- data.table::rbindlist(plotList)
  plotdat <- plotdat[order(plotdat$threshold, decreasing = TRUE), ]
   # ggplot(subset(plotdat, plotdat$subset == unique(plotdat$subset)[[4]])) +
  figure <- ggplot(plotdat) +
    geom_step(aes(x = fpr, y = sens, col = threshold)) +
    facet_wrap(~ subset, ncol = ncols) + theme_bw() +
    xlab("False Positive Rate") + ylab("True Positive Rate") +
    geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey") +
    xlim(0, 1) + ylim(0, 1) +
    scale_color_gradientn(colours = thresholdPalette, name = "Posterior") +
    ggtitle(paste("ROC Curves for Predicting", varname))
  return(figure)
}

rocTable <- function(obj, target, direction = "auto", adjust = "BH",
                     pvalue = c("wilcoxon", "logistic", "ttest", "KW"),
                     sortAUC = FALSE,...) {
  mc = match.call()
  if(!is.null(mc$minProbFilter)){
    minProbFilter = eval(mc$minProbFilter,envir=parent.frame())
  }else{
    minProbFilter = -1
  }
  notNA <- !is.na(target)
  post <- obj$posteriors[notNA, -1]
  colnames(post) <- names(obj$coefficients)

  target <- target[notNA]
  uniqVals <- unique(target)
  if(length(uniqVals) != 2 & pvalue != "KW") {
    stop("Target variable must have exactly two values!")
  }

  pvalue <- pvalue[1]

  subsets <- names(obj$posteriors)[-1]
  p <- ncol(post)
  aucs <- numeric(p)
  if(pvalue != "KW") {
    for(i in 1:p) {
      rocfit <- roc(target ~ post[, i], direction = direction)
      aucs[i] <- rocfit$auc
    }
  }

  n0 <- sum(target == uniqVals[1])
  n1 <- sum(target == uniqVals[2])
  if(pvalue == "wilcoxon") {
    pvals <- pwilcox(aucs * n0 * n1, n0, n1, lower.tail = FALSE)
  } else if(pvalue == "logistic") {
    pvals <- apply(post, 2, function(x) summary(glm(target ~ x, family = "binomial"))$coefficients[2, 4])
  } else if(pvalue == "ttest") {
    val1 <- sort(unique(target))[1]
    val2 <- sort(unique(target))[2]
    pvals <- apply(post, 2, function(x) t.test(x[target == val1], x[target == val2])$p.value)
  } else if(pvalue == "KW") {
    pvals <- apply(post, 2, function(x) kruskal.test(x, g = target)$p.value)
  } else {
    stop("unknown test type")
  }
  minpfilt = obj$levelProbs>minProbFilter
  qvals = pvals
  qvals[!minpfilt] = NA
  qvals[minpfilt] <- p.adjust(pvals[minpfilt], method = adjust)
  responseProb <- obj$levelProbs
  if(pvalue != "KW") {
    result <- data.frame(subset = subsets, responseProb  = responseProb,
                         auc = aucs, pvalue = pvals, qvalue = qvals)
  } else {
    result <- data.frame(subset = subsets, responseProb  = responseProb,
                         pvalue = pvals, qvalue = qvals)
  }

  if(sortAUC) {
    result <- result[order(result$pvalue, decreasing = TRUE), ]
  }

  return(result)
}

fdrTable <- function(obj, target) {
  notNA <- !is.na(target)
  post <- 1 - obj$posteriors[notNA, -1]
  colnames(post) <- names(obj$coefficients)
  subsets <- names(obj$coefficients)
  p <- length(subsets)

  if(is.factor(target)) {
    temp <- rep(0, length(target))
    temp[target == levels(target)[2]] <- 1
  } else if(is.numeric(target)) {
    temp <- rep(0, length(target))
    temp[target > min(target)] <- 1
  } else {
    temp <- target
  }
  target <- temp

  qvalueTable <- obj$posteriors
  plotList <- list()
  for(i in 1:p) {
    probs <- post[, i]
    postvals <- sort(unique(probs))
    estfdr <- sapply(postvals, function(x) mean(probs[probs <= x]))
    empfdr <- sapply(postvals, function(x) mean(target[probs <= x] == 0))
    power <- sapply(postvals, function(x) sum(target[probs <= x]) / sum(target))

    map <- sapply(probs, function(x) which(postvals == x))
    qvalueTable[, i + 1] <- estfdr[map]

    plotList[[i]] <- data.frame(subset = subsets[i],
                                nominal = estfdr,
                                emp = empfdr,
                                power = power)
  }

  plotList <- data.table::rbindlist(plotList)
  names(plotList)[3:4] <- c("FDR", "Power")

  results <- list()
  results$qvalues <- qvalueTable
  results$empiricalFDR <- plotList
  class(results) <- "flowReMix_fdrTable"
  return(results)
}

#' @export
#' @import ggplot2
plot.flowReMix_fdrTable <- function(x, target, subsets = NULL, varname = NULL, ncol = 5, ...) {
  plotList <- x$empiricalFDR
  if(!is.null(subsets)) {
    plotList <- subset(plotList, subset %in% subsets)
  }

  if(is.null(varname)) {
    varname <- as.character(match.call()$target)
    varname <- varname[length(varname)]
  }

  plotList <- reshape2::melt(plotList, id = c("subset", "nominal"))
  names(plotList)[3:4] <- c("Measure", "value")
  ggplot(plotList) +
    geom_line(aes(x = nominal, y = value, col = Measure, linetype = Measure)) +
    geom_abline(intercept = 0, slope = 1, col = "grey") +
    facet_wrap(~ subset, ncol = ncol) + theme_bw() +
    xlab("Nominal FDR") + ylab("Empirical") +
    ggtitle(paste("FDR and Power Curves for", varname))
}

#' @import  ggplot2
plotScatter <- function(obj, subsets = NULL,
                        target = NULL, varname = NULL,
                        ncol = 5, colPalette = NULL,
                        paletteRange = NULL,
                        summary =  "minmax") {
  if(is.null(obj$modelFrame)) {
    dat <- buildFlowFrame(obj$call, obj$data)$frame
  } else {
    dat <- obj$modelFrame
  }
  if(length(grep("%%%", dat$id[1])) == 1) {
    split <- strsplit(dat$id, "%%%")
    dat$id <- sapply(split, function(x) x[1])
    repnum <- as.numeric(sapply(split, function(x) x[2]))
    dat <- subset(dat, repnum == 1)
  }

  requireNamespace("ggplot2")
  if(!is.null(target)) {
    if(length(target) != nrow(obj$posteriors)) {
      stop("Length of target must be identical to the number of subjects in the dataset!")
    }
  }

  if(is.null(paletteRange)) {
    paletteRange <- c(0, 1)
  }

  if(is.null(colPalette)) {
    colPalette <- viridis::plasma(256)
  }

  if(is.null(varname) & !is.null(target)) {
    varname <- as.character(match.call()$target)
    varname <- varname[length(varname)]
  }

  post <- obj$posteriors
  names(post)[1] <- "id"
  post <- reshape2::melt(post, id = "id")
  names(post) <- c("id", "sub.population", "post")
  if(!is.null(target)) {
    shapes <- data.frame(id = obj$posteriors[, 1], shape = target)
    post <- merge(post, shapes, all.x = TRUE, all.y = FALSE)
  } else {
    post$shape <- NA
  }
  dat <- merge(dat, post, all.x = TRUE, all.y = FALSE,
               by.x = c("id", "sub.population"),
               by.y = c("id", "sub.population"))

  if(is.numeric(dat$treatmentvar)) {
    ctrl <- subset(dat, treatmentvar == 0)
    treat <- subset(dat, treatmentvar != 0)
  } else if(is.factor(dat$treatmentvar)) {
    ctrl <- subset(dat, treatmentvar == levels(dat$treatmentvar)[1])
    treat <- subset(dat, treatmentvar != levels(dat$treatmentvar)[1])
  } else {
    stop("treatmentvar must be numeric or factor. How did we even get this far?!")
  }

  if(summary == "minmax") {
    ctrl <- dplyr::summarise(dplyr::group_by(ctrl, sub.population, id, shape, post),
                             ctrlprop = min(prop))
    treat <- dplyr::summarise(dplyr::group_by(treat, sub.population, id, shape, post),
                              trtprop = max(prop))
  } else if(summary == "median") {
    ctrl <- dplyr::summarise(dplyr::group_by(ctrl, sub.population, id, shape, post),
                             ctrlprop = median(prop))
    treat <- dplyr::summarise(dplyr::group_by(treat, sub.population, id, shape, post),
                              trtprop = median(prop))
  } else if(summary == "mean") {
    ctrl <- dplyr::summarise(dplyr::group_by(ctrl, sub.population, id, shape, post),
                             ctrlprop = mean(prop))
    treat <- dplyr::summarise(dplyr::group_by(treat, sub.population, id, shape, post),
                              trtprop = mean(prop))
  } else {
    stop("unrecognized summary method!")
  }

  forplot <- merge(ctrl, treat)
  forplot$shape <- factor(forplot$shape)
  forplot$sub.population <- as.character(forplot$sub.population)
  if(!is.null(subsets)){
    forplot = forplot%>%filter(sub.population %in% as.character(subsets))
  }
  figure <- ggplot(forplot)
  if(is.null(target)) {
    figure <- figure + geom_point(aes(x = log(ctrlprop + 10^-9), y = log(trtprop + 10^-9),
                                      col = post))
  } else {
    figure <- figure + geom_point(aes(x = log(ctrlprop + 10^-9), y = log(trtprop + 10^-9),
                                      col = post, shape = shape)) +
      scale_shape_discrete(name = varname)
  }

  figure <- figure + scale_color_gradientn(limits = paletteRange, name = "Posterior", colors = colPalette) +
    facet_wrap(~ sub.population, scales = "free", ncol = ncol) +
    xlab("Log(Control Proportion + 10e-9) ") +
    ylab("Log(Treatment Proportion + 10e0-9)") +
    geom_abline(intercept = 0, slope = 1)

  return(figure)
}

#' @import ggplot2
plotBoxplot <- function(obj, target = NULL, varname = NULL,
                        weights = NULL, groups = c("subsets", "all"),
                        test = c("none", "logistic", "t-test", "wilcoxon", "KW"),
                        one_sided = FALSE, jitter = FALSE,
                        ncol = 5, adjust = "BH",sigdigits=5,
                        normWeights = FALSE, ...) {
  if(is.null(varname) & !is.null(target)) {
    varname <- as.character(match.call()$target)
    varname <- varname[length(varname)]
  }

  post <- obj$posteriors
  if(!is.null(target)) {
    post <- subset(post, !is.na(target))
  }

  if(is.null(weights)) {
    weights <- list()
    weights[[1]] <- rep(1, ncol(post) - 1) / (ncol(post) - 1)
    names(weights) <- c("Functionality")
    if(groups[1] == "subsets") {
      names(weights) <- c("Posterior")
    }
  } else {
    if(!is.list(weights)) {
      stop("Weights must be a named list or NULL!")
    }
    if(length(weights)>1){
      stop("Weights must be a named list of length 1")
    }
    if(is.null(names(weights))) {
      names(weights) <- paste("Score", 1:length(weights))
    }
  }

  if(groups[1] == "all") {
    groups <- list()
    groups[[1]] <- colnames(post)[-1]
    names(groups) <- c("all")
  } else if(groups[1] == "subsets") {
    groups <- lapply(colnames(post[, -1]), function(x) x)
    names(groups) <- colnames(post[-1])
  } else if(is.list(groups)) {
    if(is.null(names(groups))) {
      names(groups) <- paste("Group", 1:length(groups))
    }
  } else {
    stop("Groups must be either `all`, `subsets` or a list!")
  }

  datlist <- list()
  slot <- 1
  if(!is.null(target)) {
    target <- target[!is.na(target)]
  }

  for(i in 1:length(weights)) {
    for(j in 1:length(groups)) {
      w <- weights[[i]]
      group <- groups[[j]]
      ingroup <- which(colnames(post)[-1] %in% group)
      probs <- post[, ingroup + 1, drop = FALSE]
      w <- w[ingroup]
      if(normWeights) {
        w <- w / sum(w)
      }
      score <- apply(probs, 1, function(x) sum(x * w))

      datlist[[slot]] <- data.frame(id = post[, 1],
                                    group = names(groups)[j],
                                    measure = names(weights)[i],
                                    score = score,
                                    target = target)
      slot <- slot + 1
    }
  }
  forplot <- data.table::rbindlist(datlist)

  test <- test[1]
  if(!(test %in% c("none", "KW")) & length(unique(target)) != 2) {
    warning("Target has more than two levels, testing level 1 vs. other levels!")
  }

  facTarget <- factor(target, levels = sort(unique(target)))
  baseline <- levels(facTarget)[1]
  tempTarget <- rep(1, length(target))
  tempTarget[facTarget == baseline] <- 0
  if(test == "logistic") {
    pvalues <- by(forplot, list(forplot$group, forplot$measure),
                  function(x) summary(glm(tempTarget ~ score, family = "binomial", data = x))$coefficients[2, 4])
  } else if(test == "t-test") {
    pvalues <- by(forplot, list(forplot$group, forplot$measure),
                  function(x) t.test(x$score[tempTarget == 1], x$score[tempTarget == 0],
                                     var.equal = TRUE)$p.value)
  } else if(test == "wilcoxon") {
    pvalues <- by(forplot, list(forplot$group, forplot$measure),
                  function(x) { t <- NA
                    try(t <- wilcox.test(x$score[tempTarget == 1], x$score[tempTarget == 0],
                                     exact = FALSE)$p.value)
                    })
  } else if(test == "KW") {
    pvalues <- by(forplot, list(forplot$group, forplot$measure),
                  function(x) { t <- NA
                  try(t <- kruskal.test(x$score, target)$p.value)
                  })
  }

  if(test %in% c("t-test", "wilcoxon", "logistic", "KW")) {
    pvalues <- as.numeric(pvalues)
    if(one_sided & test != "KW") pvalues <- pvalues / 2
    adjp = p.adjust(p = pvalues, method = adjust)
    slot <- 1
    forplot$group <- as.character(forplot$group)
    for(i in 1:length(weights)) {
      for(j in 1:length(groups)) {
        index <- forplot$group == names(groups)[j] & forplot$measure == names(weights)[i]
        forplot$group[index] <- paste(forplot$group[index], ", ", test,
                                      ": ", signif(adjp[slot], sigdigits), sep = "")
        slot <- slot + 1
      }
    }
  }

  free_y <- ifelse(normWeights, "free_y", "fixed")
  figure <- ggplot(forplot)
  if(is.null(target)) {
    free_y <- ifelse(normWeights, "free_y", "fixed")
    figure <- figure + geom_boxplot(aes(x = measure, y = score),outlier.color = NA) +
      facet_wrap(~ group, scales = free_y)
    if(jitter) {
      figure <- figure + geom_point(aes(x = measure, y = score),position = position_jitterdodge(jitter.height = 0))
    }
  } else {
    figure <- figure + geom_boxplot(aes(x = measure, y = score, col = factor(target)),outlier.color = NA) +
      scale_color_discrete(name = varname)
    if(jitter) {
      figure <- figure + geom_point(aes(x = measure, y = score, col = factor(target),
                                         group = factor(target)),position = position_jitterdodge(jitter.height = 0))
    }
  }

  figure <- figure + facet_wrap(~ group, ncol = ncol, scales = free_y) +
    theme_bw() + scale_y_continuous(name = unique(forplot$measure)) +
    scale_x_discrete(name="",labels = "") + theme(axis.ticks.x = element_blank())


  return(figure)
}








