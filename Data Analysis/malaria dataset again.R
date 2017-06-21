library(flowReMix)
# Malaria dataset ----------------------------
data(malaria)
names(malaria)
table(malaria$experiment)
unique(malaria$ptid)
unique(malaria$population)
populations <- unique(malaria$population)
parents <- unique(malaria$parent)
leaves <- populations[!(populations %in% parents) ]
malaria <- subset(malaria, population %in% leaves)
unique(malaria$stim)
malaria$stimgroup[malaria$stim %in% c("PfRBC", "uRBC")] <- "RBC"
malaria$stimgroup[!(malaria$stim %in% c("PfRBC", "uRBC"))] <- "SPZ"
malaria$stim[malaria$stim == "uRBC"] <- "control"
malaria$stim <- factor(malaria$stim, levels = c("control", "PfSPZ", "PfRBC"))
isCytokine <- substring(malaria$population, nchar(malaria$population)) == "+"
malaria <- subset(malaria, isCytokine)
malaria$subset <- paste(malaria$stimgroup, "/", malaria$population, sep = "")
malaria$visitno <- factor(malaria$visitno)

malaria$infection <- TRUE
malaria$infection[malaria$ptid %in% c("60061", "50071", "20003")] <- FALSE

# Screening low counts -------------------
countlist <- by(malaria, malaria$subset, function(x) x$count)
toRemove <- sapply(countlist, function(x) mean(x > 4) < 0.05)
toRemove <- names(countlist)[toRemove]
malaria <- subset(malaria, !(subset %in% toRemove))
malaria$subset <- factor(malaria$subset)

# Analysis -----------------------
library(flowReMix)
control <- flowReMix_control(updateLag = 5, nsamp = 34, initMHcoef = 1,
                             nPosteriors = 3, centerCovariance = TRUE,
                             maxDispersion = 5000, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             lastSample = 100, isingInit = -log(89),
                             initMethod = "sparse")

tempdat <- subset(malaria, parent == "4+")
tempdat$subset <- factor(as.character(tempdat$subs))
fit <- flowReMix(cbind(count, parentcount - count) ~ visitno * stim,
                 subject_id = ptid,
                 cell_type = subset,
                 cluster_variable = visitno,
                 data = malaria,
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "sparse",
                 iterations = 10,
                 parallel = FALSE,
                 verbose = TRUE, control = control)
save(fit, file = "data analysis/results/new malaria.Robj")

# ROC analysis for infection -----------------
posteriors <- fit$posteriors
posteriors <- merge(posteriors, unique(malaria[, c(9, 14)]), all.x = TRUE)
infect <- posteriors[, ncol(posteriors)]
posteriors <- posteriors[, -c(1, ncol(posteriors))]
subsets <- colnames(posteriors)
aucs <- numeric(ncol(posteriors))
for(i in 1:ncol(posteriors)) {
  post <- posteriors[, i]
  rocfit <- roc(infect ~ post)
  aucs[i] <- rocfit$auc
}
n0 <- sum(infect)
n1 <- length(infect) - n0
pvals <- pwilcox(aucs * n0 * n1, n0, n1, lower.tail = FALSE)
qvals <- p.adjust(pvals, method = "BH")
rocResults <- data.frame(subset = subsets, auc = aucs,
                         pval = pvals, qval = qvals)

# Stability selection for graphical model --------------
assignments <- fit$assignmentList
names(assignments) <- substr(names(assignments), 1, 5)
assignments <- lapply(unique(names(assignments)), function(x) {
  do.call("rbind", assignments[names(assignments) == x])
})

reps <- 40
modelList <- list()
nsubsets <- ncol(assignments[[1]])
countCovar <- matrix(0, nrow = nsubsets, ncol = nsubsets)
for(i in 1:reps) {
  mat <- t(sapply(assignments, function(x) x[sample(1:nrow(x), 1), ]))
  colnames(mat) <- subsets
  keep <- apply(mat, 2, function(x) any(x != x[1]))
  mat <- mat[, keep]
  system.time(model <- IsingFit::IsingFit(mat, AND = TRUE, plot = TRUE))
  #coefs <- raIsing(mat, AND = TRUE, gamma = 0.9, method = "sparse")
  #plot(model)
  countCovar[keep, keep] <- countCovar[keep, keep] + (model$weiadj != 0) * sign(model$weiadj)
  #countCovar[keep, keep] <- countCovar[keep, keep] + (coefs != 0) * sign(coefs)
  print(i)
}

threshold <- 0
props <- countCovar / reps
#save(props, file = "data analysis/results/new malaria graph 4+ only.Robj")
table(props) / 2
which(props > threshold, arr.ind = TRUE)
props[abs(props) <= threshold] <- 0
sum(props != 0) / 2

# Plotting graph ---------------------
require(GGally)
library(network)
library(sna)
threshold <- 0.1
network <- props
colnames(props) <- subsets
rownames(props) <- subsets
diag(network) <- 0
# Cutting network
# network[6, 8] <- 0
# network[8, 6] <- 0
# network[6, 40] <- 0
# network[40, 6] <- 0
#######
keep <- apply(network, 1, function(x) any(abs(x) > threshold))
#keep <-  (result$qvals <= 0.05)
network <- network[keep, keep]
net <- network(props)
subsets <- names(fit$posteriors)[-1]
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
aucs <- result$auc
nodes$auc <- rocResults$auc[keep]
nodes$qvals <- rocResults$qvals[keep]

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

# Inference with connected components ----------------------
library(igraph)
network[network > 0.05] <- 1
network[network <= 0.05] <- 0
graph <- graph.adjacency(network)
comp <- components(graph)

# ROCs for vaccination
n0 <- sum(infect == 0)
n1 <- sum(infect == 1)
par(mfrow = c(1, 2))
for(i in 1:sum(comp$csize > 2)) {
  group <- which(comp$csize > 2)[i]
  group <- subsets[keep][which(comp$membership == group)]
  cols <- which(names(fit$posteriors) %in% group)
  score <- rowMeans(fit$posteriors[, cols])
  rocfit <- roc(infect ~ score)
  auc <- rocfit$auc
  pval <- pwilcox(auc * n0 * n1, n0, n1, lower.tail = FALSE)
  plot(rocfit, main = paste("Size", length(group), "AUC", round(auc, 4), "pval", round(pval, 3)))
}

