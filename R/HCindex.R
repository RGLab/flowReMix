# Helper function for computing the HC index
hcStat <- function(labels, posteriors, direction, wilcoxonPvals = NULL, a0 = 0.75) {
  if(is.null(wilcoxonPvals)) {
    wilcoxonPvals <- apply(posteriors, 2, function(x) wilcox.test(x[labels], x[!labels], exact = FALSE, alternative = "greater")$p.value)
  }
  if(direction == "greater") {
    pvals <- wilcoxonPvals
  } else if(direction == "two.sided") {
    pvals <- 2 * pmin(wilcoxonPvals, 1 - wilcoxonPvals)
  } else {
    pvals <- 1 - wilcoxonPvals
  }
  alternativeHigh <- wilcoxonPvals < 0.5
  pvals <- sort(pvals)
  pvals <- pmax(pvals, 10^-10)
  pvals <- pmin(pvals, 0.999)
  nstats <- ceiling(a0 * length(pvals))
  pvals <- pvals[1:nstats]
  HC <- (1:nstats/nstats - pvals) / sqrt(pvals * (1 - pvals))
  hcInd <- which.max(HC)
  HC <- HC[hcInd]
  index <- colnames(posteriors) %in% names(pvals)[1:hcInd]
  numLabels <- 1 - 2 * (!labels)
  index[index] <- sign(apply(posteriors[, index, drop = FALSE], 2, function(x) coef(lm(numLabels ~ x))[2]))
  index <- index / sum(abs(index))
  score <- as.numeric(as.matrix(posteriors) %*% index)
  diff <- mean(score[labels]) - mean(score[!labels])
  return(list(HC = HC, index = index, score = score, diff = diff))
}


#' A function for computing the Higher-Criticism index
#'
#' Computes the Higher-Critcism index, computes a p-value via a permutation
#' and conducts model averaging through bootstrap.
#'
#' @param fit a object of class 'flowReMix'
#'
#' @param groups the groups for which the HC index should be computed. By
#' default, the function constructs an HC index based on all cell-subsets.
#'
#' @param target a factor variable with two levels for which a classifier
#' will be constructed
#'
#' @param direction \code{character} one of two.sided, less, or greater. Direction of test.
#'
#' @param a0 Only the top a0% of cell-subsets will be considered for inclusion
#' in the HC index
#'
#' @param nPermutations \code{numeric} number of permutations for the permutation test. (default 1000)
#'
#' @param nBootstrap  \code{numeric} number of bootstrap replicates (default 200)
#'
#' @param verbose \code{logical}
#'
#' @export
hcIndex <- function(fit, groups = "all", target,
                    direction = c("two.sided", "less", "greater"), a0 = 0.75,
                    nPermutations = 1000, nBootstrap = 200,
                    verbose = TRUE) {
  # Preparing groups (if all)
  if(groups[1] == "all") {
    groups <- list()
    groups$all <- colnames(fit$posteriors)[-1]
  }

  groupSizes <- sapply(groups, length)
  if(any(groupSizes < 3)) {
    warning('Removing groups with less than 3 cell-subsets!')
    groups <- groups[groupSizes >= 3]
    if(length(groups) < 1) {
      stop("No valid groups!")
    }
  }

  # Verifying that fit is flowReMix
  if(!("flowReMix" %in% class(fit))) {
    stop("fit must be an object of class 'flowReMix'")
  }

  # Verifying that direction is valid
  direction <- direction[1]
  if(!(direction %in%  c("two.sided", "less", "greater"))) {
    stop("direction must be one of 'two.sided', 'less' or 'greater'")
  }

  # Getting target variable
  mc <- match.call()
  if(!is.null(mc$target)){
    target = mc$target
    target = enquo(target)
    subject_id = fit$subject_id
    target = fit$data %>% group_by(!!subject_id) %>% summarize(outcome=unique(!!target))%>%ungroup#%>%select(outcome)%>%unlist%>%factor
    target <- as.data.frame(target)
    post <- fit$posteriors[, 1:2]
    if(is.factor(fit$posteriors[, 1])) {
      target[, 1] <- factor(target[, 1], levels = levels(fit$posteriors[, 1]))
    }
    target <- merge(post, target, sort = FALSE)
    target <- target[, 3]
  } else {
    stop("target variable must be specified!")
  }

  # Validating target
  if(!is.factor(target)) {
    if(length(unique(target)) != 2) {
      stop("target must have exactly two-values")
    } else {
      target <- factor(target, levels = unique(target))
      warning("target is not a factor, converting to factor, verify that test direction is correct.")
    }
  } else {
    if(length(levels(target)) != 2) {
      stop("target must have exactly two levels!")
    }
  }

  # preparing labels
  baseline <- levels(target)[1]
  alternative <- levels(target)[2]
  labels <- target == alternative

  # preparing posterior probability table
  posteriors <- aggregate$posteriors
  nGroups <- length(groups)
  inGroups <- unlist(groups) %>% unique()
  inGroups <- colnames(posteriors) %in% inGroups
  posteriors <- posteriors[, inGroups]
  groupIndices <- lapply(groups, function(g) colnames(posteriors) %in% g)
  anyAllSame <- apply(posteriors, 2, function(x) all(x == x[1]))
  if(any(anyAllSame)) {
    stop("Some cell-subsets have constant posterior probabilities,
         remove those from the groups and try again.")
  }

  # computing observed HC index
  pvals <- apply(posteriors, 2, function(x) wilcox.test(x[labels], x[!labels], exact = FALSE, alternative = "greater")$p.value)
  observed <- lapply(groupIndices, function(g) {
    hcStat(labels, posteriors[, g], direction, wilcoxonPvals = pvals[g], a0 = a0)
  })

  # permutation test
  M <- nPermutations
  differences <- matrix(nrow = M, ncol = length(groups))
  hcStats <- matrix(nrow = M, ncol = length(groups))
  if(verbose) {
    print("Conducting permutation test")
    pb <- txtProgressBar(min = 0, max = M, style = 3)
  }
  for(i in 1:M) {
    permLabels <- labels[order(runif(length(labels)))]
    pvals <- apply(posteriors, 2, function(x) wilcox.test(x[permLabels], x[!permLabels], exact = FALSE, alternative = "greater")$p.value)
    permFit <- lapply(groupIndices, function(g) {
      hcStat(permLabels, posteriors[, g], direction, wilcoxonPvals = pvals[g], a0 = a0)
    })
    differences[i, ] <- sapply(permFit, function(x) x$diff)
    hcStats[i, ] <- sapply(permFit, function(x) x$HC)
    if(verbose) setTxtProgressBar(pb, i)
  }
  if(verbose) close(pb)

  # computing pvalues based on permutations
  indexPvals <- numeric(nGroups)
  hcPvals <- numeric(nGroups)
  obsDiffs <- sapply(observed, function(x) x$diff)
  obsHC <- sapply(observed, function(x) x$HC)
  for(i in 1:nGroups) {
    indexPvals[i] <- mean(differences[, i] >= obsDiffs[i])
    hcPvals[i] <- mean(hcStats[, i] >= obsHC[i])
  }
  names(indexPvals) <- names(groups)
  names(hcPvals) <- names(groups)

  # Bootstrapping
  M <- nBootstrap
  inclusions <- matrix(0, ncol = length(groups), nrow = ncol(posteriors))
  rownames(inclusions) <- colnames(posteriors)
  colnames(inclusions) <- names(groups)
  nBaseline <- sum(target == baseline)
  nAlt <- sum(target == alternative)
  if(verbose) {
    print("Model Averaging")
    pb <- txtProgressBar(min = 0, max = M, style = 3)
  }
  for(i in 1:M) {
    baselineSamp <- sample(which(target == baseline), nBaseline, TRUE)
    alternativeSamp <- sample(which(target == alternative), nAlt, TRUE)
    bootSamp <- c(baselineSamp, alternativeSamp)
    bootlabels <- labels[bootSamp]
    bootpost <- posteriors[bootSamp, ]
    pvals <- apply(bootpost, 2, function(x) wilcox.test(x[bootlabels], x[!bootlabels], exact = FALSE, alternative = "greater")$p.value)
    bootFit <- lapply(groupIndices, function(g) {
      hcStat(bootlabels, bootpost[, g], direction, wilcoxonPvals = pvals[g], a0 = a0)
    })
    for(j in 1:nGroups) {
      inclusions[groupIndices[[j]], j] <- inclusions[groupIndices[[j]], j] + sign(bootFit[[j]]$index)
    }
    if(verbose) setTxtProgressBar(pb, i)
  }
  if(verbose) close(pb)

  result <- list(indices = observed,
                 indexPvals = indexPvals,
                 hcPvals = hcPvals,
                 propensities = inclusions / M)
  class(result) <- "flowReMix_HC"
  return(result)
}
