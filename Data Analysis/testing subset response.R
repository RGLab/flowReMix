library(flowReMix)
library(pROC)
cummean <- function(x) cumsum(x) / 1:length(x)
data(rv144)
#set.seed(502)
omit <- paste("P", c(1001, 1013, 1019, 1023, 1031, 1034, 1039, 1045,
                     1060, 1095, 1099, 1100, 1109, 1177, 1180, 1187,
                     1201, 1215, 1216, 1224, 1227, 1232, 1242, 1284),
              sep = "")
par(mfrow = c(1, 1), mar = rep(4, 4))
data <- rv144
data <- subset(data, !(ptid %in% omit))
leaves <- unique(data$population)
selected_populations = c(c(1, 2, 3, 5, 4, 6, 7))
selected_populations = c(c(1, 2, 5, 3, 6, 7))
#selected_populations = c(c(1, 2, 5))
data <- subset(data, population %in% leaves[selected_populations])
data$population <- factor(data$population)
data <- subset(data, stim != "sebctrl")
data$treatment <- as.numeric(data$stim == "env")
data$ptid <- as.numeric(data$ptid)
data$ptid[data$vaccine == "VACCINE"] <- data$ptid[data$vaccine == "VACCINE"] * 10^4
data$prop <- data$count / data$parentcount
data$population <- as.factor(data$populataion)
data <- data[order(data$population, data$ptid, data$stim, decreasing = FALSE), ]

preAssign <- function(dat) {
  subsets <- unique(dat$population)
  nSubsets <- length(subsets)
  preAssign <- numeric(nSubsets)
  prop <- dat$count / dat$parentcount
  for(j in 1:nSubsets) {
    negctrl <- prop[dat$stim == "negctrl" & dat$population == subsets[j]]
    env <- prop[dat$stim == "env" & dat$population == subsets[j]]
    preAssign[j] <- ifelse(env > negctrl, -1, 0)
  }
  result <- data.frame(id = dat$ptid[1], subset = subsets, assign = preAssign)
  return(result)
}
preAssignment <- by(data, INDICES = data$ptid, preAssign)
preAssignment <- do.call("rbind", preAssignment)

vaccine <- as.numeric(by(data, data$ptid, function(x) x$vaccine[1] == "VACCINE"))
system.time(fit <- subsetResponseMixtureRcpp(count ~  treatment,
                                         sub.population = factor(data$population),
                                         N = parentcount, id =  ptid, treatment = treatment,
                                         data = data, preAssignment = preAssignment,
                                         randomAssignProb = 0.0,
                                         weights = NULL,
                                         rate = 1, updateLag = 5, nsamp = 30, maxIter = 25,
                                         sparseGraph = TRUE, betaDiserpsion = TRUE,
                                         covarianceMethod = c("dense"),
                                         centerCovariance = FALSE,
                                         initMHcoef = 3))
#save(fit, file = "dispersed model 2.Robj")
#save(fit, file = "results/binom model.Robj")
#save(fit, file = "results/dispersed model 2 wAssignment.Robj")
#load("results/dispersed model 2 wAssignment.Robj")

require(pROC)
vaccine <- as.vector(by(data, INDICES = data$ptid, FUN = function(x) x$vaccine[1] == "VACCINE"))
posteriors <- fit$posteriors[, 2:ncol(fit$posteriors), drop = FALSE]
par(mfrow = c(2, 3), mar = rep(3, 4))
for(i in 1:length(selected_populations)) {
  rocfit <- roc(vaccine ~ posteriors[, i])
  print(plot(rocfit, main = paste(leaves[selected_populations[i]], "- AUC", round(rocfit$auc, 3))))
}

par(mfrow = c(2, 3), mar = rep(3, 4))
for(i in 1:length(selected_populations)) {
  post <- posteriors[, i]
  treatment <- vaccine[order(post)]
  uniquePost <- sort(unique(post))
  nominalFDR <- sapply(uniquePost, function(x) mean(post[post <= x]))
  empFDR <- sapply(uniquePost, function(x) 1 - mean(vaccine[post <= x]))
  power <- sapply(uniquePost, function(x) sum(vaccine[post <= x]) / sum(vaccine))
  print(plot(nominalFDR, empFDR, type = "l", xlim = c(0, 1), ylim = c(0, 1), col = "red", main = leaves[selected_populations[i]]))
  lines(nominalFDR, power, col = "blue", lty = 2)
  abline(a = 0, b = 1)
  abline(v = c(0.05, 0.1), h = c(0.8, 0.9), col = "grey")
}

posteriors <- fit$posteriors[, 2:ncol(fit$posteriors), drop = FALSE]
forplot <- list()
vaccine <- as.vector(by(data, INDICES = data$ptid, FUN = function(x) x$vaccine[1] == "VACCINE"))
for(i in 1:length(selected_populations)) {
  post <- posteriors[, i]
  negprop <- log(data$count / data$parentcount)[data$population == leaves[selected_populations[i]] & data$stim == "negctrl"]
  envprop <- log(data$count / data$parentcount)[data$population == leaves[selected_populations[i]] & data$stim == "env"]
  forplot[[i]] <- data.frame(subset = leaves[selected_populations[i]],
                             negprop = negprop, envprop = envprop,
                             posterior = 1 - post, vaccine = vaccine)
}

forplot <- do.call("rbind", forplot)
require(ggplot2)
print(ggplot(forplot) +
  geom_point(aes(x = negprop, y = envprop, col = posterior, shape = vaccine)) +
  facet_wrap(~ subset, scales = 'free') +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + scale_colour_gradientn(colours=rainbow(4)))



# Booleans dataset --------------------------------------------
data("rv144_booleans")
bySubset <- by(data.frame(booleans$stim, booleans$nonstim), booleans$Subset, function(x) x)
largerThanThershold <- sapply(bySubset, function(x) colSums(x >5))

require(reshape2)
booldata <- melt(booleans, c("PTID", "Subset"))
names(booldata)[3:4] <- c("stim", "count")

forParentcount <- rv144
forParentcount <- as.data.frame(forParentcount)
forParentcount <- subset(forParentcount,
                         parent == "4+" & stim == "env")
forParentcount <- forParentcount[, c(2, 5, 11)]
forParentcount <- unique(forParentcount)

booldata <- merge(booldata, forParentcount, by.x = "PTID", by.y = "ptid",
                  all.x = TRUE, all.y = FALSE)
booldata <- subset(booldata, Subset != "!TNFa&!IFNg&!IL4&!IL2&!CD154&!IL17a")
booldata$treatment <- as.numeric(booldata$stim == "stim")
uniquepop <- unique(booldata$Subset)
booldata <- with(booldata, booldata[order(Subset, PTID, stim, decreasing = FALSE), ])
booldata <- subset(booldata, !is.na(Subset))
allsubset <- booldata
booldata <- with(booldata, booldata[order(Subset, PTID, stim, decreasing = FALSE), ])

vaccine <- as.numeric(by(booldata, booldata$PTID, function(x) x$vaccine[1] == "VACCINE"))
require(pROC)
system.time(fit <- subsetResponseMixtureRcpp(count ~  treatment,
                                             sub.population = factor(booldata$Subset),
                                             N = parentcount, id =  PTID,
                                             data = booldata,
                                             treatment = treatment,
                                             randomAssignProb = 0,
                                             weights = NULL,
                                             rate = 1, updateLag = 10,
                                             nsamp = 100, maxIter = 30,
                                             initMHcoef = 3,
                                             covarianceMethod = "sparse",
                                             sparseGraph = TRUE,
                                             centerCovariance = FALSE))
#save(fit, file = "results/boolean dispersed fit.Robj")
#save(fit, file = "results/boolean dispersed fit2.Robj")
load("results/boolean dispersed fit2.Robj")
subsetIndex <- 1:length(subsets)
#subsetIndex <- c(3, 5, 6, 10, 12, 13, 14, 20, 23)
subsets <- unique(booldata$Subset)

require(pROC)
posteriors <- fit$posteriors[, 2:ncol(fit$posteriors), drop = FALSE]
#par(mfrow = c(4, 6), mar = rep(1, 4))
par(mfrow = c(3, 3), mar = rep(1, 4))
auc <- numeric(length(subsets))
for(j in 1:length((subsetIndex))) {
  i <- subsetIndex[j]
  try(rocfit <- roc(!vaccine ~ posteriors[, i]))
  try(rocfit <- roc(vaccine ~ posteriors[, i]))
  auc[i] <- rocfit$auc
  print(plot(rocfit, main = paste(i, "- AUC", round(rocfit$auc, 3))))
}


par(mfrow = c(4, 6), mar = rep(2, 4))
par(mfrow = c(3, 3), mar = rep(2, 4))
for(j in 1:length((subsetIndex))) {
  i <- subsetIndex[j]
  post <- posteriors[, i]
  treatment <- vaccine[order(post)]
  uniquePost <- sort(unique(post))
  nominalFDR <- sapply(uniquePost, function(x) mean(post[post <= x]))
  empFDR <- sapply(uniquePost, function(x) 1 - mean(vaccine[post <= x]))
  power <- sapply(uniquePost, function(x) sum(vaccine[post <= x]) / sum(vaccine))
  print(plot(nominalFDR, empFDR, type = "l", xlim = c(0, 1), ylim = c(0, 1), col = "red",
             main = paste(i, "- AUC ", round(auc[i], 2))))
  lines(nominalFDR, power, col = "blue", lty = 2)
  abline(a = 0, b = 1)
  abline(v = c(0.05, 0.1), h = c(0.8, 0.9), col = "grey")
}


forplot <- list()
for(j in 1:length(subsetIndex)) {
  i <- subsetIndex[j]
  post <- posteriors[, i]
  negprop <- log(booldata$count / booldata$parentcount)[booldata$Subset == subsets[i] & booldata$stim == "nonstim"]
  envprop <- log(booldata$count / booldata$parentcount)[booldata$Subset == subsets[i] & booldata$stim == "stim"]
  forplot[[i]] <- data.frame(subset = subsets[i],
                             negprop = negprop, envprop = envprop,
                             posterior = 1 - post, vaccine = vaccine)
}

forplot <- do.call("rbind", forplot[c(1:21, 23)])
require(ggplot2)
print(ggplot(forplot) +
        geom_point(aes(x = negprop, y = envprop, col = posterior, shape = vaccine == 1)) +
        facet_wrap(~ subset, scales = 'free', ncol = 6) +
        geom_abline(slope = 1, intercept = 0) +
        theme_bw() + scale_colour_gradientn(colours=rainbow(4)))



library(gridExtra)
table <- data.frame(index = 1:length(subsets), subset = subsets,
                    responseProb = round(fit$levelProbs,2),
                    AUC = round(auc, 2))
pdf("figures/cell subset table B.pdf", height=8, width=6)
grid.table(table, rows = NULL)
dev.off()

par(mfrow = c(1, 1))
plot(fit$isingfit)

# Plotting network -----------------------
require(GGally)
library(network)
library(sna)
network <- fit$isingfit$weiadj
nodes <- ggnet2(network(net))$data
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
edges$width <- with(edges, abs(width) / max(abs(width)) * sign(width))
nodes$probs <- fit$levelProbs
coefficients <- sapply(fit$coefficients, function(x) x[2])
for(i in 1:nrow(nodes)) {
  if(coefficients[i] < 0) {
    nodes$probs[i] <- 1 - nodes$probs[i]
  }
}
nodes$auc <- auc
nodes$mFunctional <- 1:nrow(nodes) >= 6
names(edges)[5] <- "Dependence"
ggplot() +
  scale_colour_gradient(limits=c(0, 1), low="white", high = "black") +
  geom_segment(data = edges, aes(x = xstart, y = ystart,
                                          xend = xend, yend = yend,
                                          alpha = Dependence,
                                          col = Dependence),
               size = 1) +
  scale_fill_gradientn(colours=rainbow(4)) +
  geom_point(data = nodes[1:5, ], aes(x = x, y = y, fill = auc), shape = 25,
             size = 8, col = "white") +
  geom_point(data = nodes[6:23, ], aes(x = x, y = y, fill = auc), shape = 21,
             size = 8, col = "white") +
  scale_size(range = c(0.3, 1)) +
  scale_shape(solid = FALSE) +
  geom_text(data = nodes, aes(x = x, y = y, label = 1:23)) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())





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
malaria <- subset(malaria, stim %in% c("control", "PfRBC"))
malaria <- subset(malaria, population %in% leaves)

# just one subset
parents <- parents[parents %in% unique(malaria$parent)]
fitList <- list()
for(j in 11:length(parents)) {
  # Fitting model --------------------------
  i <- j
  tempdat <- subset(malaria, parent == parents[j])
  tempdat <- tempdat[order(tempdat$ptid, tempdat$population, tempdat$stim), ]
  by(tempdat, tempdat$population, function(x) do.call("cbind", by(x, x$visitno, function(y) y$count)))
  if(i == 1) {
    tempdat <- subset(tempdat, population != "4+/IL4+")
  }
  if(i == 3) {
    tempdat <- subset(tempdat, !(population %in% c("PD-1+/IL21+", "PD-1+/IL4+")))
  }
  #else if(i == 4) {
  #   next
  # } else if(i == 6) {
  #   tempdat <- subset(tempdat, population != "8+/CXCR5+/TNFa+")
  #   next
  # } else if(i == 7) {
  #   next
  # } else if(j == 9) {
  #   next
  # } else if(j == 10) {
  #   next
  # }
  tempdat$treatment <- factor(as.numeric(tempdat$stim == "PfRBC"))
  system.time(fit <- subsetResponseMixtureRcpp(count ~  treatment * visitno,
                                               sub.population = factor(tempdat$population),
                                               N = parentcount, id =  ptid,
                                               data = tempdat,
                                               treatment = treatment,
                                               updateLag = 5,
                                               nsamp = 40, maxIter = 15,
                                               initMHcoef = 1,
                                               randomAssignProb = 0.1,
                                               covarianceMethod = "sparse",
                                               sparseGraph = TRUE, betaDispersion = TRUE,
                                               centerCovariance = FALSE))
  fitList[[j]] <- fit

  # Scatter and counts over time plots
  tempdat <- subset(malaria, parent == parents[j])
  tempdat$daynum <- 0
  tempdat$daynum[tempdat$visitno == "Day 9"] <- 1
  tempdat$daynum[tempdat$visitno == "pos"] <- 2
  tempdat$daynum[tempdat$visitno == "Day 28"] <- 3
  tempdat$daynum[tempdat$visitno == "Day 56"] <- 4
  tempdat$daynum[tempdat$visitno == "Day 168"] <- 5
  tempdat$logprop <- with(tempdat, log((count + 1)/parentcount + 10^-5))
  tempdat$visitno <- factor(tempdat$visitno,
                            levels = c("Day 0", "Day 9", "pos", "Day 28", "Day 56", "Day 168"))

  ggplot(tempdat) +
    geom_boxplot(aes(x = visitno, y = logprop, col = stim)) +
    theme_bw() +
    facet_wrap(~ population, scales = "free")

  control <- subset(tempdat, stim == "control")
  stim <- subset(tempdat, stim == "PfRBC")
  names(stim)[13] <- "stimprop"
  names(control)[13] <- "controlprop"
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
    pred <- data.frame(times = rep(times, 2),
                       type = c(rep("nonresponse", ntimes),
                                rep("response", ntimes)))
    coefs <- fit$coefficients[[i]]
    # nonresponse
    predictions <- rep(coefs[[1]], ntimes)
    predictions <- predictions + c(0, coefs[c(6, 4, 5, 3)])
    pred$eta[1:ntimes] <- predictions
    # response
    predictions <- predictions + coefs[c(2, 10, 8, 9, 7)]
    pred$eta[(ntimes + 1):(ntimes * 2)] <- predictions
    pred$prop <- 1 / (1 + exp(-pred$eta))
    pred$population <- names(fit$coefficients)[i]
    predlist[[i]] <- pred
  }

  predicted <- do.call("rbind", predlist)
  ggplot(subset(predicted, population != "PD-1+/IFNg+")) +
    geom_line(aes(x = as.numeric(times), y = log(prop), col = population, linetype = type))

  # Jackknife errors --------------------
  replicateDataset <- function(data, replicate) {
    data$ptid <- paste(data$ptid, "%%%", replicate, sep = "")
    return(data)
  }

  jackFitList <- list()
  for(i in 1:length(unique(tempdat$ptid))) {
    id <- unique(tempdat$ptid)[i]
    jackdat <- subset(tempdat, ptid != id)
    jackdat$visitno <- factor(as.character(jackdat$visitno))
    jackdat$treatment <- (as.numeric(jackdat$stim == "PfRBC"))
    #jackdat <- do.call("rbind", lapply(1:7, function(x) replicateDataset(jackdat, x)))
    jackfit <- subsetResponseMixtureRcpp(count ~  treatment * visitno,
                                         sub.population = factor(jackdat$population),
                                         N = parentcount, id =  ptid,
                                         data = jackdat,
                                         treatment = treatment,
                                         updateLag = 3,
                                         nsamp = 20, maxIter = 2,
                                         initMHcoef = 1,
                                         randomAssignProb = 0.1,
                                         covarianceMethod = "sparse",
                                         sparseGraph = TRUE, betaDispersion = FALSE,
                                         centerCovariance = FALSE,
                                         dataReplicates = 2)
    jackFitList[[i]] <- jackfit
  }
}

par(mfrow = c(2, 3), mar = rep(3, 4))
for(i in 1:length(fitList)) {
  if(!is.null(fitList[[i]])){
    plot(fitList[[i]]$isingfit)
  }
}


# Plotting graphs
par(mfrow = c(2, 3), mar = rep(3, 4))
for(i in 1:length(fitList)) {
  if(!is.null(fitList[[i]])){
    plot(fitList[[i]]$isingfit)
  }
}

# Creating dataset for plotting pairwise scatters.
j <- 5
fit <- fitList[[j]]
tempdat <- subset(malaria, parent == parents[j])
tempdat$daynum <- 0
tempdat$daynum[tempdat$visitno == "Day 9"] <- 1
tempdat$daynum[tempdat$visitno == "pos"] <- 2
tempdat$daynum[tempdat$visitno == "Day 28"] <- 3
tempdat$daynum[tempdat$visitno == "Day 56"] <- 4
tempdat$daynum[tempdat$visitno == "Day 168"] <- 5
tempdat$logprop <- with(tempdat, log((count + 1)/parentcount + 10^-5))
control <- subset(tempdat, stim == "control")
stim <- subset(tempdat, stim == "PfRBC")
names(stim)[13] <- "stimprop"
names(control)[13] <- "controlprop"
tempdat <- merge(stim, control, by = c("ptid", "population", "visitno"))

ggplot(tempdat) + geom_point(aes(x = controlprop, y = stimprop, col = ptid)) +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(population ~ visitno, scales = "free_y")

ggplot(tempdat, aes(x = daynum.x, y = stimprop - controlprop, col = ptid)) +
  geom_point() + geom_line() +
  facet_wrap(~ population, scales = "free_y") +
  geom_hline(yintercept = 0)

# Extract predicted values
times <- unique(tempdat$visitno)
ntimes <- length(times)
predlist <- list()
for(i in 1:length(fit)) {
  pred <- data.frame(times = rep(times, 2),
                     type = c(rep("response", ntimes),
                              rep("nonresponse", ntimes)))
  coefs <- fit$coefficients[[i]]
  # nonresponse
  predictions <- rep(coefs[[1]], ntimes)
  predictions <- predictions + c(0, coefs[c(6, 4, 5, 3)])
  pred$eta[1:ntimes] <- predictions
  # response
  predictions <- predictions + coefs[c(2, 10, 8, 9, 7)]
  pred$eta[(ntimes + 1):(ntimes * 2)] <- predictions
  pred$prop <- 1 / (1 + exp(-pred$eta))
  pred$population <- names(fit$coefficients)[i]
  predlist[[i]] <- pred
}

predicted <- do.call("rbind", predlist)
ggplot(subset(predicted, population != "PD-1+/IFNg+")) +
  geom_line(aes(x = as.numeric(times), y = log(prop), col = population, linetype = type))






