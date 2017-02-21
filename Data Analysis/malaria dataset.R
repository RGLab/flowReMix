preAssign <- function(dat) {
  subsets <- unique(dat$population)
  nSubsets <- length(subsets)
  if(dat$infection[1]) {
    preAssign <- rep(-1, nSubsets)
  } else {
    preAssign <- rep(0, nSubsets)
  }
  result <- data.frame(id = dat$ptid[1], subset = subsets, assign = preAssign)
  return(result)
}

# Malaria dataset ----------------------------
data(malaria)
names(malaria)
table(malaria$experiment)
unique(malaria$ptid)
unique(malaria$population)
populations <- unique(malaria$population)
parents <- unique(malaria$parent)
leaves <- populations[!(populations %in% parents) ]
unique(malaria$stim)
malaria <- subset(malaria, stim %in% c("uRBC", "PfRBC"))
malaria <- subset(malaria, population %in% leaves)
malaria$infection <- TRUE
malaria$infection[malaria$ptid %in% c("60061", "50071", "20003")] <- FALSE

# just one subset
parents <- parents[parents %in% unique(malaria$parent)]
parents <- parents[c(1, 2, 3, 5, 6, 9, 10, 11)]
fitList <- list()
pvalList <- list()
analysisorder <- c(7, 1, 2, 3, 4, 5, 6, 8)
for(j in 1:length(parents)) {
  # Fitting model -------------------------
  i <- j
  tempdat <- subset(malaria, parent == parents[j])
  tempdat <- tempdat[order(tempdat$ptid, tempdat$population, tempdat$stim), ]
  tempdat$treatment <- factor(as.numeric((tempdat$stim == "PfRBC") * tempdat$infection))

  countList <- by(tempdat, tempdat$population, function(x) do.call("cbind", by(x, x$visitno, function(y) y$count)))
  removeSubset <- sapply(countList, function(x) mean(x >= 2) < 0.05)
  toRemove <- names(removeSubset)[removeSubset]
  tempdat <- subset(tempdat, !(tempdat$population %in% toRemove))

  preAssignment <- by(tempdat, INDICES = tempdat$ptid, preAssign)
  preAssignment <- do.call("rbind", preAssignment)

  system.time(fit <- subsetResponseMixtureRcpp(count ~  treatment * visitno,
                                               sub.population = factor(tempdat$population),
                                               N = parentcount, id =  ptid,
                                               data = tempdat,
                                               treatment = treatment,
                                               updateLag = 5,
                                               nsamp = 30, maxIter = 12,
                                               initMHcoef = 1,
                                               randomAssignProb = 0.0,
                                               covarianceMethod = "sparse",
                                               smallCounts = TRUE,
                                               sparseGraph = TRUE,
                                               betaDispersion = FALSE,
                                               centerCovariance = FALSE,
                                               preAssignment = preAssignment,
                                               dataReplicates = 3))
  fitList[[j]] <- fit
  #save(fit, file = paste("data analysis/results/", sep = "", gsub("/", "", x = parents[j]), "malaria binom glmnet.Robj"))
  load(file = paste("data analysis/results/", sep = "", gsub("/", "", x = parents[j]), "malaria binom glmnet.Robj"))
  #load(file = "data analysis/results/4+ malaria binom.Robj")
  #load(file = "data analysis/results/4+CXCR5+ malaria binom.Robj")

  # Scatter and counts over time plots
  tempdat <- subset(malaria, parent == parents[j])
  tempdat <- subset(tempdat, !(tempdat$population %in% toRemove))
  tempdat <- tempdat[order(tempdat$ptid, tempdat$population, tempdat$stim), ]
  tempdat$treatment <- factor(as.numeric((tempdat$stim == "PfRBC") * tempdat$infection))
  # if(j == 3) {
  #   tempdat <- subset(tempdat, !(population %in% c("PD-1+/IL4+", "PD-1+/IL21+")))
  # }
  tempdat$daynum <- 0
  tempdat$daynum[tempdat$visitno == "Day 9"] <- 1
  tempdat$daynum[tempdat$visitno == "pos"] <- 2
  tempdat$daynum[tempdat$visitno == "Day 28"] <- 3
  tempdat$daynum[tempdat$visitno == "Day 56"] <- 4
  tempdat$daynum[tempdat$visitno == "Day 168"] <- 5
  tempdat$logprop <- with(tempdat, log(count/parentcount + 10^-4))
  tempdat$visitno <- factor(tempdat$visitno,
                            levels = c("Day 0", "Day 9", "pos", "Day 28", "Day 56", "Day 168"))

  ggplot(tempdat) +
    geom_boxplot(aes(x = visitno, y = logprop, col = stim)) +
    theme_bw() +
    facet_wrap(~ population, scales = "free")

  control <- subset(tempdat, stim == "uRBC")
  stim <- subset(tempdat, stim == "PfRBC")
  names(stim)[15] <- "stimprop"
  names(control)[15] <- "controlprop"
  tempdat <- merge(stim, control, by = c("ptid", "population", "visitno"))

  ggplot(tempdat) + geom_point(aes(x = controlprop, y = stimprop, col = ptid)) +
    geom_abline(intercept = 0, slope = 1) +
    facet_grid(visitno ~ population, scales = "free")

  ggplot(tempdat) +
    geom_line(alpha = 0.5, aes(x = daynum.x, y = stimprop - controlprop, col = ptid)) +
    facet_wrap(~ population, scales = "free_y") +
    geom_hline(yintercept = 0)

  ggplot(tempdat) +
    geom_boxplot(alpha = 0.5, aes(x = visitno, y = stimprop - controlprop)) +
    facet_wrap(~ population, scales = "free_y") +
    geom_hline(yintercept = 0)

  # Predicted values ----------------------
  times <- unique(tempdat$visitno)
  ntimes <- length(times)
  predlist <- list()
  for(i in 1:length(fit$coefficients)) {
    pred <- data.frame(times = sort(times))
    coefs <- fit$coefficients[[i]]
    # nonresponse
    predictions <- rep(coefs[[1]], ntimes)
    predictions <- predictions + c(0, coefs[c(6, 7, 4, 5, 3)])
    pred$nulleta <- predictions
    # response
    predictions <- predictions + rep(coefs[2], ntimes)
    predictions <- predictions + c(0, coefs[c(12, 11, 9, 10, 8)])
    pred$alteta <- predictions
    pred$nullprop <- 1 / (1 + exp(-pred$nulleta))
    pred$altprop <- 1 / (1 + exp(-pred$alteta))
    pred$population <- names(fit$coefficients)[i]
    predlist[[i]] <- pred
  }

  predicted <- do.call("rbind", predlist)

  ggplot(tempdat) +
    facet_wrap( ~ population, scales = "free_y", nrow = 3) +
    geom_line(data = tempdat, alpha = 0.35, aes(x = daynum.x, y = stimprop, linetype = ptid), col = "red") +
    geom_line(data = predicted, aes(x = as.numeric(times) - 1, y = log(altprop + 10^-4)), col = "red") +
    geom_line(data = tempdat, alpha = 0.35, aes(x = daynum.x, y = controlprop, linetype = ptid), col = "blue") +
    geom_line(data = predicted, aes(x = as.numeric(times) - 1, y = log(nullprop + 10^-4)), col = "blue") +
    theme_bw()

  # Jackknife errors --------------------
  jackFitList <- list()
  tempdat <- subset(malaria, parent == parents[j])
  tempdat <- subset(tempdat, !(tempdat$population %in% toRemove))
  tempdat <- tempdat[order(tempdat$ptid, tempdat$population, tempdat$stim), ]
  tempdat$treatment <- factor(as.numeric((tempdat$stim == "PfRBC") * tempdat$infection))
  for(i in 1:length(unique(tempdat$ptid))) {
    id <- unique(tempdat$ptid)[i]
    jackdat <- subset(tempdat, ptid != id)
    jackdat$visitno <- factor(as.character(jackdat$visitno))

    preAssignment <- by(jackdat, INDICES = jackdat$ptid, preAssign)
    preAssignment <- do.call("rbind", preAssignment)

    jackfit <- subsetResponseMixtureRcpp(count ~  treatment * visitno,
                                         sub.population = factor(jackdat$population),
                                         N = parentcount, id =  ptid,
                                         data = jackdat,
                                         treatment = treatment,
                                         updateLag = 5,
                                         nsamp = 30, maxIter = 12,
                                         initMHcoef = 1,
                                         randomAssignProb = 0.1,
                                         preAssignment = preAssignment,
                                         covarianceMethod = "sparse",
                                         sparseGraph = TRUE, betaDispersion = FALSE,
                                         centerCovariance = FALSE,
                                         dataReplicates = 5)
    jackFitList[[i]] <- jackfit
  }
  #save(jackFitList, file = paste("data analysis/results/", sep = "", gsub("/", "", x = parents[j]), "malaria binom glmnet jackknife.Robj"))

  # Jackknife estimate of curve-difference variability
  nsubsets <- length(jackFitList[[1]]$coefficients)
  times <- unique(tempdat$visitno)
  ntimes <- length(times)
  zvals <- numeric(nsubsets)
  zvalMat <- matrix(nrow = nsubsets, ncol = ntimes)
  for(i in 1:nsubsets) {
    reps <- length(jackFitList)
    diffMat <- matrix(nrow = reps, ncol = ntimes)
    for(j in 1:reps) {
      jackfit <- jackFitList[[j]]
      coefs <- jackfit$coefficients[[i]]
      # nonresponse
      nonresponse <- rep(coefs[[1]], ntimes)
      nonresponse <- nonresponse + c(0, coefs[c(6, 7, 4, 5, 3)])
      # response
      response <- nonresponse + rep(coefs[[2]], ntimes)
      response <- response + c(0, coefs[c(11, 12, 9, 10, 8)])
      diffMat[j, ] <- response - nonresponse
    }
    meanDiff <- colMeans(diffMat)
    varMat <- t(apply(diffMat, 1, function(x) (x - meanDiff)))
    varMat <- t(varMat) %*% varMat * (reps - 1) / reps
    ivec <- rep(1, ntimes) / ntimes
    zvals[i] <- mean(meanDiff) / as.numeric(sqrt(t(ivec) %*% varMat %*% ivec))
    zvalMat[i, ] <- meanDiff / sqrt(diag(varMat))
  }
  pvals <- pt(zvals, df = 11, lower.tail = FALSE)
  names(pvals) <- names(jackfit$coefficients)
  pvalMat <- pt(zvalMat, df = 11, lower.tail = FALSE)
  bhMat <- matrix(p.adjust(pvalMat, method = "BH"), ncol = ncol(pvalMat), byrow = FALSE)
  pvalList[[j]] <- pvals
  p.adjust(pvals, method = "BH")
}


# Results for parent subset analysis ----------------------------------
populations <- unique(malaria$population)
parents <- unique(malaria$parent)
leaves <- populations[!(populations %in% parents) ]
parents <- parents[parents %in% unique(malaria$parent)]
parents <- parents[c(1, 2, 3, 5, 6, 9, 10, 11)]
times <- unique(malaria$visitno)
ntimes <- length(times)
pvalList <- list()
for(j in 1:length(parents)) {
  load(file = paste("data analysis/results/", sep = "", gsub("/", "", x = parents[j]), "malaria binom glmnet jackknife.Robj"))
  nsubsets <- length(jackFitList[[1]]$coefficients)
  zvals <- numeric(nsubsets)
  zvalMat <- matrix(nrow = nsubsets, ncol = ntimes)
  for(i in 1:nsubsets) {
    reps <- length(jackFitList)
    diffMat <- matrix(nrow = reps, ncol = ntimes)
    for(k in 1:reps) {
      jackfit <- jackFitList[[k]]
      coefs <- jackfit$coefficients[[i]]
      # nonresponse
      nonresponse <- rep(coefs[[1]], ntimes)
      nonresponse <- nonresponse + c(0, coefs[c(6, 7, 4, 5, 3)])
      # response
      response <- nonresponse + rep(coefs[[2]], ntimes)
      response <- response + c(0, coefs[c(11, 12, 9, 10, 8)])
      diffMat[k, ] <- response - nonresponse
    }
    meanDiff <- colMeans(diffMat)
    varMat <- t(apply(diffMat, 1, function(x) (x - meanDiff)))
    varMat <- t(varMat) %*% varMat * (reps - 1) / reps
    ivec <- rep(1, ntimes) / ntimes
    zvals[i] <- mean(meanDiff) / as.numeric(sqrt(t(ivec) %*% varMat %*% ivec))
    zvalMat[i, ] <- meanDiff / sqrt(diag(varMat))
  }
  pvals <- pt(zvals, df = 11, lower.tail = FALSE)
  names(pvals) <- names(jackfit$coefficients)
  pvalMat <- pt(zvalMat, df = 11, lower.tail = FALSE)
  bhMat <- matrix(p.adjust(pvalMat, method = "BH"), ncol = ncol(pvalMat), byrow = FALSE)
  pvalList[[j]] <- pvals
  p.adjust(pvals, method = "BH")
}

names(pvalList) <- parents
level <- 0.05
qvalList <- lapply(pvalList, function(x) p.adjust(x, method = "BH"))
global <- sapply(qvalList, function(x) any(x < level))
qvalList <- qvalList[global]
subset <- lapply(qvalList, function(x) x < sum(global) / length(global) * level)





# Joint analysis -----------------------------
parents <- unique(malaria$parent)
leaves <- populations[!(populations %in% parents) ]
parents <- parents[parents %in% unique(malaria$parent)]
parents <- parents[c(1, 2, 3, 5, 6, 9, 10, 11)]
tempdat <- subset(malaria, parent %in% parents)
tempdat <- tempdat[order(tempdat$ptid, tempdat$population, tempdat$stim), ]
tempdat$treatment <- factor(as.numeric((tempdat$stim == "PfRBC") * tempdat$infection))
countList <- by(tempdat, tempdat$population, function(x) do.call("cbind", by(x, x$visitno, function(y) y$count)))
removeSubset <- sapply(countList, function(x) mean(x >= 2) < 0.05)
toRemove <- names(removeSubset)[removeSubset]
tempdat <- subset(tempdat, !(tempdat$population %in% toRemove))

preAssignment <- by(tempdat, INDICES = tempdat$ptid, preAssign)
preAssignment <- do.call("rbind", preAssignment)

system.time(fit <- subsetResponseMixtureRcpp(count ~  treatment * visitno,
                                             sub.population = factor(tempdat$population),
                                             N = parentcount, id =  ptid,
                                             data = tempdat,
                                             treatment = treatment,
                                             updateLag = 10,
                                             nsamp = 40, maxIter = 30,
                                             initMHcoef = 1,
                                             randomAssignProb = 1,
                                             covarianceMethod = "sparse",
                                             smallCounts = TRUE,
                                             sparseGraph = TRUE,
                                             betaDispersion = FALSE,
                                             centerCovariance = FALSE,
                                             preAssignment = preAssignment,
                                             dataReplicates = 10))
save(fit, "Data Analysis/results/malaria everything.Robj")

