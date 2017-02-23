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
treatmentStim <- "PfRBC"
controlStim <- "uRBC"
malaria <- subset(malaria, stim %in% c(controlStim, treatmentStim))
malaria <- subset(malaria, population %in% leaves)
malaria$infection <- TRUE
malaria$infection[malaria$ptid %in% c("60061", "50071", "20003")] <- FALSE

# just one subset ------------------
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
  tempdat$treatment <- factor(as.numeric((tempdat$stim == treatmentStim) * tempdat$infection))

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

  control <- subset(tempdat, stim == controlStim)
  stim <- subset(tempdat, stim == treatmentStim)
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

tempdat$visit <- tempdat$visitno
tempdat$visitno[!tempdat$infection] <- "Day 0"
tempdat$treatment <- as.numeric(tempdat$stim == treatmentStim)

countList <- by(tempdat, tempdat$population, function(x) do.call("cbind", by(x, x$visitno, function(y) y$count)))
removeSubset <- sapply(countList, function(x) mean(x >= 2) < 0.05)
toRemove <- names(removeSubset)[removeSubset]
tempdat <- subset(tempdat, !(tempdat$population %in% toRemove))
leaves <- unique(tempdat$population)
cytokines <- sapply(as.character(leaves), function(x) stringr::str_sub(x, start = -1, end = -1) == "+")
cytokines <- leaves[cytokines]
tempdat <- subset(tempdat, population %in% cytokines)

system.time(fit <- subsetResponseMixtureRcpp(count ~  treatment * visitno,
                                             sub.population = factor(tempdat$population),
                                             N = parentcount, id =  ptid,
                                             data = tempdat,
                                             treatment = treatment,
                                             updateLag = 10,
                                             nsamp = 50, maxIter = 30,
                                             initMHcoef = 1,
                                             randomAssignProb = 0,
                                             covarianceMethod = "sparse",
                                             regressionMethod = "sparse",
                                             isingMethod = "sparse",
                                             centerCovariance = FALSE,
                                             preAssignment = NULL,
                                             dataReplicates = 5))
#save(fit, file = "Data Analysis/results/malaria everything cytokines 4.Robj")
#load(file = "Data Analysis/results/malaria everything cytokines.Robj")
#plot(fit$isingfit)

# Better estimation of graphical model ----------------------
assignments <- fit$assignmentList
names(assignments) <- sapply(names(assignments), stringr::str_sub, 1, 5)
assignments <- lapply(unique(names(assignments)),
                      function(x) do.call("rbind", assignments[names(assignments) == x]))
reps <- 1000
modelList <- list()
nsubsets <- ncol(assignments[[1]])
countCovar <- matrix(0, nrow = nsubsets, ncol = nsubsets)
for(i in 1:reps) {
  mat <- t(sapply(assignments, function(x) x[sample(1:nrow(x), 1), ]))
  keep <- apply(mat, 2, function(x) any(x != x[1]))
  mat <- mat[, keep]
  model <- IsingFit::IsingFit(mat, AND = FALSE, plot = FALSE)
  modelList[[i]] <- model
  #plot(model)
  countCovar[keep, keep] <- countCovar[keep, keep] + (model$weiadj != 0) * sign(model$weiadj)
  print(i)
}

props <- countCovar / reps
table(props)
threshold <- 0.07
which(props > 0.1, arr.ind = TRUE)
props[abs(props) <= threshold] <- 0
sum(props != 0) / 2

# Plotting graph ---------------------
require(GGally)
library(network)
library(sna)
network <- props
keep <- apply(network, 1, function(x) any(abs(x) > 0))
network <- network[keep, keep]
net <- network(props)
nodes <- ggnet2(network, label = cytokines[keep])$data
edges <- matrix(nrow = sum(network != 0)/2, ncol = 5)
p <- nrow(network)
row <- 1
for(j in 2:p) {
  for(i in 1:(j-1)) {
    if(network[i, j] != 0) {
      edges[row, ] <- unlist(c(nodes[i, 6:7], nodes[j, 6:7], network[i, j]))
      row <- row + 1
    }
  }
}

edges <- data.frame(edges)
names(edges) <- c("xstart", "ystart", "xend", "yend", "width")
coefficients <- sapply(fit$coefficients, function(x) x[2])
nodes$probs <- fit$levelProbs[keep]

names(edges)[5] <- "Dependence"
lims <- max(abs(props))
ggplot() +
  scale_colour_gradient2(limits=c(-lims, lims), low="dark red", high = "dark green") +
  geom_segment(data = edges, aes(x = xstart, y = ystart,
                                 xend = xend, yend = yend,
                                 col = Dependence,
                                 alpha = Dependence),
               size = 1) +
  scale_fill_gradient2(low = "white", high = "red") +
  #scale_fill_gradientn(colours = rainbow(4))+
  geom_point(data = nodes, aes(x = x, y = y, fill = probs), shape = 21,
             size = 12, col = "grey") +
  scale_size(range = c(0.3, 1)) +
  scale_shape(solid = FALSE) +
  geom_text(data = nodes, aes(x = x, y = y, label = nodes$label), size = 1.8) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(), legend.position = "none")


# Prediction plots -------------------------
parents <- sapply(cytokines, function(x) {
  split <- strsplit(x, "/")[[1]]
  return(paste0(split[1:(length(split) - 1)], sep = "", collapse = "/"))
})
parents <- unlist(parents)
cytType <- unlist(sapply(cytokines, function(x) {
  split <- strsplit(x, "/")[[1]]
  return(paste(split[length(split)], sep = ""))
}))

ind <- 1
plotdat <- tempdat
plotdat <- subset(plotdat, population %in% cytokines[parents == unique(parents)[ind]])
#plotdat <- subset(plotdat, population %in% sigSubsets)
plotdat$daynum <- 0
plotdat$daynum[plotdat$visit == "Day 9"] <- 1
plotdat$daynum[plotdat$visit == "pos"] <- 2
plotdat$daynum[plotdat$visit == "Day 28"] <- 3
plotdat$daynum[plotdat$visit == "Day 56"] <- 4
plotdat$daynum[plotdat$visit == "Day 168"] <- 5
plotdat$logprop <- with(plotdat, log(count/parentcount + 10^-4))
plotdat$visit <- factor(plotdat$visit,
                          levels = c("Day 0", "Day 9", "pos", "Day 28", "Day 56", "Day 168"))

ggplot(plotdat) +
  geom_boxplot(aes(x = visit, y = logprop, col = stim)) +
  theme_bw() +
  facet_wrap(~ population, scales = "free")

control <- subset(plotdat, stim == "uRBC")
stim <- subset(plotdat, stim == "PfRBC")
names(stim)[16] <- "stimprop"
names(control)[16] <- "controlprop"
plotdat <- merge(stim, control, by = c("ptid", "population", "visit"))

# Predicted values ----------------------
coefList <- fit$coefficients[names(fit$coefficients) %in% unique(plotdat$population)]
times <- unique(plotdat$visit)
ntimes <- length(times)
predlist <- list()
for(i in 1:length(coefList)) {
  pred <- data.frame(times = sort(times))
  coefs <- coefList[[i]]
  # nonresponse
  predictions <- rep(coefs[[1]], ntimes)
  predictions <- predictions + c(0, coefs[c(6, 7, 4, 5, 3)])
  pred$nulleta <- predictions
  # response
  predictions <- predictions + rep(coefs[2], ntimes)
  predictions <- predictions + c(0, coefs[c(11, 12, 9, 10, 8)])
  pred$alteta <- predictions
  pred$nullprop <- 1 / (1 + exp(-pred$nulleta))
  pred$altprop <- 1 / (1 + exp(-pred$alteta))
  pred$population <- names(coefList)[i]
  predlist[[i]] <- pred
}

predicted <- do.call("rbind", predlist)

require(ggplot2)
ggplot(plotdat) +
  facet_wrap( ~ population, scales = "free_y", nrow = 3) +
  geom_line(data = plotdat, alpha = 0.35, aes(x = daynum.x, y = stimprop, group = ptid, linetype = infection.x), col = "red") +
  geom_line(data = predicted, aes(x = as.numeric(times) - 1, y = log(altprop + 10^-4)), col = "red") +
  geom_line(data = plotdat, alpha = 0.35, aes(x = daynum.x, y = controlprop, group = ptid, linetype = infection.x), col = "blue") +
  geom_line(data = predicted, aes(x = as.numeric(times) - 1, y = log(nullprop + 10^-4)), col = "blue") +
  theme_bw()



# Jackknife inference ----------------------
jackFitList <- list()
for(i in 1:length(unique(tempdat$ptid))) {
  id <- unique(tempdat$ptid)[i]
  jackdat <- subset(tempdat, ptid != id)
  jackdat$visitno <- factor(as.character(jackdat$visitno))

  system.time(jackfit <- subsetResponseMixtureRcpp(count ~  treatment * visitno,
                                               sub.population = factor(jackdat$population),
                                               N = parentcount, id =  ptid,
                                               data = jackdat,
                                               treatment = treatment,
                                               updateLag = 10,
                                               nsamp = 50, maxIter = 20,
                                               initMHcoef = 1,
                                               randomAssignProb = 1,
                                               covarianceMethod = "sparse",
                                               regressionMethod = "sparse",
                                               isingMethod = "sparse",
                                               centerCovariance = FALSE,
                                               preAssignment = NULL,
                                               dataReplicates = 5))
  jackFitList[[i]] <- jackfit
}
#save(jackFitList, file = "Data Analysis/results/malaria everything cytokines jackknife 2.Robj")
#load(file = "Data Analysis/results/malaria everything cytokines jackknife.Robj")

nsubsets <- length(jackFitList[[1]]$coefficients)
zvals <- numeric(nsubsets)
beZvals <- numeric(nsubsets)
bePvals <- numeric(nsubsets)
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


  beDiff <- as.numeric(diffMat %*% c(0 , 0, -1, 1, 0, 0))
  meanBE <- mean(beDiff)
  varBE <- sum((beDiff - meanBE)^2) * (reps - 1) / reps
  beZvals[i] <- meanBE / sqrt(varBE)
  bePvals[i] <- pt(beZvals[i], lower.tail = FALSE, df = reps - 1)
}
pvals <- pt(zvals, df = 11, lower.tail = FALSE)
names(pvals) <- names(jackfit$coefficients)
qvals <- p.adjust(pvals, method = "BH")

# Preparing qvalue table
parents <- sapply(cytokines, function(x) {
  split <- strsplit(x, "/")[[1]]
  return(paste0(split[1:(length(split) - 1)], sep = "", collapse = "/"))
})
parents <- unlist(parents)
cytType <- unlist(sapply(cytokines, function(x) {
  split <- strsplit(x, "/")[[1]]
  return(paste(split[length(split)], sep = ""))
}))

qvalMat <- matrix(nrow = length(unique(cytType)), ncol = length(unique(parents)))
qvalMat <- data.frame(qvalMat)
names(qvalMat) <- unique(parents)
rownames(qvalMat) <- unique(cytType)
for(i in 1:length(qvals)) {
  row <- which(rownames(qvalMat) == cytType[i])
  col <- which(names(qvalMat) == parents[i])
  qvalMat[row, col] <- qvals[i]
}

table <- round(qvalMat, 3)
for(i in 1:nrow(table)) {
  for(j in 1:ncol(table)) {
    if(is.na(table[i, j])) next
    if(table[i, j] < 0.1) {
      table[i, j] <- paste0("\\colorbox{yellow}{", table[i, j], "}")
    }
  }
}

print(xtable::xtable(table), sanitize.text.function = function(x) x)




