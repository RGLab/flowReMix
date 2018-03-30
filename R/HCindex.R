# Helper function for computing the HC index
hcStat <- function(labels, posteriors, direction, wilcoxonPvals = NULL, a0 = 0.75, plus = TRUE) {
  if(is.null(wilcoxonPvals)) {
    wilcoxonPvals <- apply(posteriors, 2, function(x) wilcox.test(x[labels], x[!labels], exact = FALSE, alternative = "greater")$p.value)
  }
  N <- length(wilcoxonPvals)
  if(direction == "greater") {
    pvals <- wilcoxonPvals
  } else if(direction == "two.sided") {
    pvals <- 2 * pmin(wilcoxonPvals, 1 - wilcoxonPvals)
  } else {
    pvals <- 1 - wilcoxonPvals
  }
  alternativeHigh <- wilcoxonPvals < 0.5
  pvalOrder <- order(pvals)
  pvals <- sort(pvals)
  # pvals <- pmax(pvals, 10^-10)
  # pvals <- pmin(pvals, 0.999)
  nstats <- ceiling(a0 * N)
  pvals <- pvals[1:nstats]
  unifQuants <- 1:nstats/N
  HC <- (unifQuants - pvals) / sqrt(pvals * (1 - pvals))
  if(plus) {
    lbound <- -1
  } else {
    lbound <- 1 / N
  }
  if(any(pvals >= lbound)) {
    hcInd <- which.max(HC[pvals >= lbound]) + sum(pvals < lbound)
  } else {
    hcInd <- length(HC)
  }

  HC <- HC[hcInd]
  # index <- colnames(posteriors) %in% names(pvals)[1:hcInd]
  index <- 1:N %in% pvalOrder[1:hcInd]
  numLabels <- 1 - 2 * (!labels)
  index[index] <- sign(apply(posteriors[, index, drop = FALSE], 2, function(x) coef(lm(numLabels ~ x))[2]))
  index[is.na(index)] <- 0
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
#' @param a0 Only the top a0% of cell-subsets will be considered for inclusion
#' in the HC index
#'
#' @export
hcIndex <- function(fit, groups = "all", target,
                    direction = c("two.sided", "less", "greater"), a0 = 0.75,
                    nPermutations = 1000, nBootstrap = 200,
                    plus = TRUE, innovated = FALSE, scale = FALSE,
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

  # if innovated then can only do one group
  if(innovated & length(groups) > 1) {
    stop("innovated HC can only be computed for a single group at a time!")
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

  # removing NA targets
  posteriors <- aggregate$posteriors
  posteriors <- posteriors[!is.na(target), ]
  target <- target[!is.na(target)]
  if(is.factor(target)) {
    levels(target) <- levels(target)[!is.na(levels(target))]
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
  nGroups <- length(groups)
  inGroups <- unlist(groups) %>% unique()
  inGroups <- colnames(posteriors) %in% inGroups
  posteriors <- posteriors[, inGroups]
  if(!innovated) {
    groupIndices <- lapply(groups, function(g) colnames(posteriors) %in% g)
  } else {
    groupIndices <- list()
    groupIndices$iHC <- 1:min(nrow(posteriors), ncol(posteriors))
  }
  anyAllSame <- apply(posteriors, 2, function(x) all(x == x[1]))
  if(any(anyAllSame)) {
    stop("Some cell-subsets have constant posterior probabilities,
         remove those from the groups and try again.")
  }

  # computing observed HC index
  if(innovated) {
    tpost <- scale(posteriors, scale = scale)
    eigen <- svd(tpost)
    tpost <- eigen$u %*% diag(eigen$d)
  } else {
    tpost <- posteriors
  }
  pvals <- apply(tpost, 2, function(x) wilcox.test(x[labels], x[!labels], exact = FALSE, alternative = "greater")$p.value)
  observed <- lapply(groupIndices, function(g) {
    hcStat(labels, tpost[, g], direction, wilcoxonPvals = pvals[g], a0 = a0,
           plus = plus)
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
    pvals <- apply(tpost, 2, function(x) wilcox.test(x[permLabels], x[!permLabels], exact = FALSE, alternative = "greater")$p.value)
    permFit <- lapply(groupIndices, function(g) {
      hcStat(permLabels, tpost[, g], direction, wilcoxonPvals = pvals[g],
             a0 = a0, plus = plus)
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
    if(innovated) {
      tpost <- scale(bootpost, scale = scale)
      eigen <- svd(tpost)
      tpost <- eigen$u %*% diag(eigen$d)
    } else {
      tpost <- posteriors
    }
    pvals <- apply(tpost, 2, function(x) wilcox.test(x[bootlabels], x[!bootlabels], exact = FALSE, alternative = "greater")$p.value)
    bootFit <- lapply(groupIndices, function(g) {
      hcStat(bootlabels, tpost[, g], direction, wilcoxonPvals = pvals[g],
             a0 = a0, plus = plus)
    })
    if(!innovated) {
      for(j in 1:nGroups) {
        inclusions[groupIndices[[j]], j] <- inclusions[groupIndices[[j]], j] + sign(bootFit[[j]]$index)
      }
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
