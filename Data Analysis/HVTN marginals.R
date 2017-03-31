# Loading Data --------------------------------
# hvtn <- read.csv(file = "data/merged_505_stats.csv")
# names(hvtn) <- tolower(names(hvtn))
# hvtn <- subset(hvtn, !is.na(ptid))
# saveRDS(hvtn, file = "data/505_stats.rds")

# Getting marginals -----------------------------
library(flowReMix)
hvtn <- readRDS(file = "data/505_stats.rds")
length(unique(hvtn$name))
length(unique(hvtn$ptid))
length(unique(hvtn$population))
unique(hvtn$population)
unique(hvtn$stim)
nchars <- nchar(as.character(unique(hvtn$population)))
#marginals <- unique(hvtn$population)[nchars < 26]
marginals <- unique(hvtn$population)[nchars == 26]
marginals <- subset(hvtn, population %in% marginals)
marginals <- subset(marginals, stim %in% c("negctrl", "CMV", "VRC ENV A",
                                           "VRC ENV B", "VRC ENV C",
                                           "VRC GAG B", "VRC NEF B",
                                           "VRC POL 1 B", "VRC POL 2 B"))
marginals <- subset(marginals, !(population %in% c("4+", "8+")))
marginals <- subset(marginals, !(population %in% c("8+/107a-154-IFNg-IL2-TNFa-", "4+/107a-154-IFNg-IL2-TNFa-")))
marginals$stim <- factor(as.character(marginals$stim))
marginals$population <- factor(as.character(marginals$population))

# Descriptives -------------------------------------
library(ggplot2)
marginals$prop <- marginals$count / marginals$parentcount
# ggplot(marginals) + geom_boxplot(aes(x = population, y = log(prop), col = stim))

require(dplyr)
negctrl <- subset(marginals, stim == "negctrl")
negctrl <- summarize(group_by(negctrl, ptid, population), negprop = mean(prop))
negctrl <- as.data.frame(negctrl)
marginals <- merge(marginals, negctrl, all.x = TRUE)

# ggplot(subset(marginals, stim != "negctrl" & parent == "4+")) +
#   geom_point(aes(x = log(negprop), y = log(prop)), size = 0.25) +
#   facet_grid(stim ~ population, scales = "free") +
#   theme_bw() +
#   geom_abline(intercept = 0, slope = 1)

# Setting up data for analysis ---------------------------
unique(marginals$stim)
gag <- subset(marginals, stim %in% c("VRC GAG B", "negctrl"))
gag$subset <- factor(paste("gag", gag$population, sep = "/"))
gag$stimGroup <- "gag"
pol <-subset(marginals, stim %in% c("negctrl", "VRC POL 1 B", "VRC POL 2 B"))
pol$subset <- factor(paste("pol", pol$population, sep = "/"))
pol$stimGroup <- "pol"
cmv <- subset(marginals, stim %in% c("negctrl", "CMV"))
cmv$subset <- factor(paste("cmv", cmv$population, sep = "/"))
cmv$stimGroup <- "cmv"
env <- subset(marginals, stim %in% c("negctrl", "VRC ENV C", "VRC ENV B", "VRC ENV A"))
env$subset <- factor(paste("env", env$population, sep = "/"))
env$stimGroup <- "env"
nef <- subset(marginals, stim %in% c("negctrl", "VRC NEF B"))
nef$subset <- factor(paste("nef", nef$population, sep = "/"))
nef$stimGroup <- "nef"
subsetDat <- rbind(gag, pol, cmv, env, nef)
subsetDat$stim <- as.character(subsetDat$stim)
subsetDat$stim[subsetDat$stim == "negctrl"] <- 0
subsetDat$stim <- factor(subsetDat$stim)

# Getting outcomes -------------------------------
treatmentdat <- read.csv(file = "data/rx_v2.csv")
names(treatmentdat) <- tolower(names(treatmentdat))
treatmentdat$ptid <- factor(gsub("-", "", (treatmentdat$ptid)))
treatmentdat <- subset(treatmentdat, ptid %in% unique(subsetDat$ptid))

# Fitting the model ------------------------------
library(flowReMix)
control <- flowReMix_control(updateLag = 15, nsamp = 250, initMHcoef = 1,
                             nPosteriors = 2, centerCovariance = TRUE,
                             maxDispersion = 10^3 / 2, minDispersion = 10^8,
                             randomAssignProb = 0.2, intSampSize = 50,
                             initMethod = "binom")

subsetDat$batch <- factor(subsetDat$batch..)
subsetDat$stimGroup <- factor(subsetDat$stimGroup)
fit <- flowReMix(cbind(count, parentcount - count) ~ stim,
                 subject_id = ptid,
                 cell_type = subset,
                 cluster_variable = stim,
                 data = subsetDat,
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "betabinom",
                 iterations = 20,
                 verbose = TRUE, control = control)

#save(fit, file = "Data Analysis/results/HVTN betabinom 2.Robj")
#load(file = "Data Analysis/results/HVTN binom stim response.Robj")
#load(file = "Data Analysis/results/HVTN binom w batch sparse 2.Robj")
#load(file = "Data Analysis/results/HVTN binom 1.Robj")
#load(file = "Data Analysis/results/HVTN betabinom 1.Robj")


# ROC plots -----------------------------
require(pROC)
outcome <- treatmentdat[, c(13, 15)]
posteriors <- fit$posteriors
posteriors <- merge(posteriors, outcome,
                    by.x = "ptid", by.y = "ptid",
                    all.x = TRUE)
ctrlCol <- ncol(posteriors)
par(mfrow = c(4, 5), mar = rep(3, 4))
aucs <- numeric(ctrlCol - 2)
n1 <- sum(posteriors$control == 1)
n0 <- sum(posteriors$control == 0)
for(i in 2:(ctrlCol - 1)) {
  post <- posteriors[, i]
  outcome <- posteriors[, ctrlCol]
  rocfit <- roc(outcome ~ post)
  subset <- names(posteriors)[i]
  auc <- rocfit$auc
  try(plot(rocfit, main = paste(subset, round(auc, 3))))
  aucs[i - 1] <- auc
}

pvals <- pwilcox(aucs * n0 * n1, n0, n1, lower.tail = FALSE)
pvals <- 2 * pmin(pvals,  1 - pvals)
qvals <- p.adjust(pvals, method = "BY")
subsets <- names(posteriors)[-c(1, ctrlCol)]
result <- data.frame(subsets, aucs, pvals, qvals)
result[order(result$aucs, decreasing = TRUE), ]

# Scatter plots -----------------------
require(reshape2)
outcome <- treatmentdat[, c(13, 15)]
posteriors <- fit$posteriors
posteriors <- melt(posteriors, id.vars = c("ptid"))
names(posteriors)[2:3] <- c("subset", "posterior")
posteriors <- merge(posteriors, outcome,
                    by.x = "ptid", by.y = "ptid",
                    all.x = TRUE)
forplot <- merge(subsetDat, posteriors, all.x = TRUE,
                 by.x = c("ptid", "subset"),
                 by.y = c("ptid", "subset"))

ggplot(subset(forplot, stim != "0" & parent == "4+")) +
  geom_point(aes(x = log(negprop), y = log(prop), col = 1 - posterior), size = 0.25) +
  facet_grid(stim ~ population, scales = "free") +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_gradientn(colours=rainbow(4))


# Stability selection for graphical model ------------------------
assignments <- fit$assignmentList
names(assignments) <- substr(names(assignments), 1, 9)
assignments <- lapply(unique(names(assignments)), function(x) {
  do.call("rbind", assignments[names(assignments) == x])
})
reps <- 20
modelList <- list()
nsubsets <- ncol(assignments[[1]])
countCovar <- matrix(0, nrow = nsubsets, ncol = nsubsets)
for(i in 1:reps) {
  mat <- t(sapply(assignments, function(x) x[sample(1:nrow(x), 1), ]))
  colnames(mat) <- names(fit$coefficients)
  keep <- apply(mat, 2, function(x) any(x != x[1]))
  mat <- mat[, keep]
  model <- IsingFit::IsingFit(mat, AND = TRUE, plot = TRUE)
  modelList[[i]] <- model
  #plot(model)
  countCovar[keep, keep] <- countCovar[keep, keep] + (model$weiadj != 0) * sign(model$weiadj)
  print(i)
}

props <- countCovar / reps
table(props)
threshold <- 0.5
which(props > threshold, arr.ind = TRUE)
props[abs(props) <= threshold] <- 0
sum(props != 0) / 2
#save(props, file = "Data Analysis/results/HVTN betabinom 2 graph.Robj")


# Plotting graph ---------------------
require(GGally)
library(network)
library(sna)
network <- props
keep <- apply(network, 1, function(x) any(abs(x) >= threshold)) | (aucs >= 0.7)
network <- network[keep, keep]
net <- network(props)
subsets <- names(fit$coefficients)
nodes <- ggnet2(network, label = subsets[keep])$data
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
nodes$auc <- aucs[keep]
nodes$qvals <- result$qvals[keep]
nodes$sig <- nodes$qvals < 0.1

names(edges)[5] <- "Dependence"
lims <- max(abs(props))
library(ggplot2)
ggplot() +
  scale_colour_gradient2(limits=c(-lims, lims), low="dark red", high = "dark green") +
  geom_segment(data = edges, aes(x = xstart, y = ystart,
                                 xend = xend, yend = yend,
                                 col = Dependence,
                                 alpha = abs(Dependence)),
               size = 1) +
  #scale_fill_gradient2(low = "white", high = "red", limits = c(0.7, 1)) +
  scale_fill_gradientn(colours = rainbow(4))+
  geom_point(data = nodes, aes(x = x, y = y, fill = auc), shape = 21,
             size = 8, col = "grey") +
  scale_shape(solid = FALSE) +
  geom_text(data = nodes, aes(x = x, y = y, label = nodes$label), size = 1.8) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

nodes[order(nodes$auc),]

# Grouping by stim ---------------------------
assignments <- fit$assignmentList
names(assignments) <- substr(names(assignments), 1, 9)
assignments <- lapply(unique(names(assignments)), function(x) {
  do.call("rbind", assignments[names(assignments) == x])
})
sapply(assignments, function(samp) mean(apply(samp, 1, function(x) sum(x == 0) < 40)))
subsets <- names(fit$coefficients)
summarizeBySubset <- function(samp, subsets, threshold = 0) {
  #colnames(samp) <- substring(subsets, 5) # for cell-type groups
  colnames(samp) <- substring(subsets, 1, 6) # for stim groups
  result <- sapply(unique(colnames(samp)), function(x) {
    index <- colnames(samp) == x
    t(apply(samp, 1, function(y) sum(y[index]) > threshold))
  })
}
collapsed <- lapply(assignments, summarizeBySubset, subsets, threshold = 0)

# AUC
posteriors <- matrix(nrow = length(collapsed), ncol = ncol(collapsed[[1]]))
for(i in 1:length(collapsed)) posteriors[i, ] <- colMeans(collapsed[[i]])
posteriors <- data.frame(posteriors)
posteriors$ptid <- fit$posteriors[, 1]
colnames(posteriors) <- c(colnames(collapsed[[1]]), "ptid")

outcome <- treatmentdat[, c(13, 15)]
posteriors <- merge(posteriors, outcome,
                    by.x = "ptid", by.y = "ptid",
                    all.x = TRUE)
ctrlCol <- ncol(posteriors)
par(mfrow = c(4, 5), mar = rep(3, 4))
aucs <- numeric(ctrlCol - 2)
n1 <- sum(posteriors$control == 1)
n0 <- sum(posteriors$control == 0)
for(i in 2:(ctrlCol - 1)) {
  post <- 1 - posteriors[, i]
  outcome <- posteriors[, ctrlCol]
  rocfit <- roc(outcome ~ post)
  subset <- names(posteriors)[i]
  auc <- rocfit$auc
  try(plot(rocfit, main = paste(subset, round(auc, 3))))
  aucs[i - 1] <- auc
}

pvals <- pwilcox(aucs * n0 * n1, n0, n1, lower.tail = FALSE)
pvals <- 2 * pmin(pvals,  1 - pvals)
qvals <- p.adjust(pvals, method = "BY")
subsets <- names(posteriors)[-c(1, ctrlCol)]
result <- data.frame(subsets, aucs, pvals, qvals)
result[order(result$aucs, decreasing = TRUE), ]


# Graph
cells <- colnames(collapsed[[1]])
assignments <- collapsed
reps <- 40
modelList <- list()
nsubsets <- ncol(assignments[[1]])
countCovar <- matrix(0, nrow = nsubsets, ncol = nsubsets)
for(i in 1:reps) {
  mat <- t(sapply(assignments, function(x) x[sample(1:nrow(x), 1), ]))
  colnames(mat) <- names(cells)
  keep <- apply(mat, 2, function(x) any(x != x[1]))
  mat <- mat[, keep]
  model <- IsingFit::IsingFit(mat, AND = TRUE, plot = TRUE)
  modelList[[i]] <- model
  #plot(model)
  countCovar[keep, keep] <- countCovar[keep, keep] + (model$weiadj != 0) * sign(model$weiadj)
  print(i)
}

props <- countCovar / reps
#save(props, file = "Data Analysis/results/HVTN binom w batch sparse 2 graph.Robj")

table(props)
threshold <- 0.1
which(props > threshold, arr.ind = TRUE)
props[abs(props) <= threshold] <- 0
sum(props != 0) / 2

# Plotting graph
require(GGally)
library(network)
library(sna)
network <- props
keep <- apply(network, 1, function(x) any(abs(x) >= threshold)) | TRUE
network <- network[keep, keep]
net <- network(props)
subsets <- names(fit$coefficients)
nodes <- ggnet2(network, label = cells)$data
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
nodes$auc <- aucs[keep]
nodes$qvals <- result$qvals[keep]
nodes$sig <- nodes$qvals < 0.1

names(edges)[5] <- "Dependence"
lims <- max(abs(props))
library(ggplot2)
ggplot() +
  scale_colour_gradient2(limits=c(-lims, lims), low="dark red", high = "dark green") +
  geom_segment(data = edges, aes(x = xstart, y = ystart,
                                 xend = xend, yend = yend,
                                 col = Dependence,
                                 alpha = abs(Dependence)),
               size = 1) +
  #scale_fill_gradient2(low = "white", high = "red", limits = c(0.7, 1)) +
  scale_fill_gradientn(colours = rainbow(4))+
  geom_point(data = nodes, aes(x = x, y = y, fill = auc), shape = 21,
             size = 8, col = "grey") +
  scale_shape(solid = FALSE) +
  geom_text(data = nodes, aes(x = x, y = y, label = nodes$label), size = 1.8) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

