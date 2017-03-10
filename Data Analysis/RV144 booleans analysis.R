require(pROC)
require(reshape2)
data("rv144_booleans")
bySubset <- by(data.frame(booleans$stim, booleans$nonstim), booleans$Subset, function(x) x)
largerThanThershold <- sapply(bySubset, function(x) colSums(x >5))

booldata <- melt(booleans, c("PTID", "Subset"))
names(booldata)[3:4] <- c("stim", "count")

booldata <- by(booldata, INDICES = list(booldata$PTID, booldata$stim), function(x) {
  x$parentcount <- sum(x$count)
  return(x)
})
booldata <- do.call("rbind", booldata)

booldata <- subset(booldata, Subset != "!TNFa&!IFNg&!IL4&!IL2&!CD154&!IL17a")
booldata$treatment <- as.numeric(booldata$stim == "stim")
uniquepop <- unique(booldata$Subset)
booldata <- with(booldata, booldata[order(Subset, PTID, stim, decreasing = FALSE), ])
booldata <- subset(booldata, !is.na(Subset))
allsubset <- booldata
booldata <- with(booldata, booldata[order(Subset, PTID, stim, decreasing = FALSE), ])

# Naming ------------------
subsets <- unique(booldata$Subset)
booldata$Subset <- as.character(booldata$Subset)
nfunctions <- numeric(length(subsets))
for(i in 1:length(subsets)) {
  split <- strsplit(as.character(subsets[i]), "&")[[1]]
  first <- substr(split, 1, 1)
  nfunction <- sum(first != "!")
  nfunctions[i] <- nfunction
  name <- paste(split[first != "!"], collapse = ",")
  booldata$nfunction[booldata$Subset == subsets[[i]]] <- nfunction
  booldata$Subset[booldata$Subset == subsets[[i]]] <- name
}
subsets <- unique(booldata$Subset)
booldata <- with(booldata, booldata[order(Subset, PTID, stim, decreasing = FALSE), ])

# Getting vaccine information --------------------
data("rv144")
rv144 <- rv144[order(rv144$ptid), ]
vaccine <- as.numeric(by(rv144, rv144$ptid, function(x) x$vaccine[1] == "VACCINE"))
vaccine <- vaccine[unique(rv144$ptid) %in% unique(booldata$PTID)]

# Analysis -------------
system.time(fit <- subsetResponseMixtureRcpp(count ~  treatment,
                                             sub.population = factor(booldata$Subset),
                                             N = parentcount, id =  PTID,
                                             data = booldata,
                                             treatment = treatment,
                                             randomAssignProb = 0,
                                             weights = NULL,
                                             updateLag = 20, nsamp = 50, maxIter = 30,
                                             initMHcoef = 3,
                                             covarianceMethod = "sparse",
                                             isingMethod = "sparse",
                                             regressionMethod = "betabinom",
                                             centerCovariance = FALSE,
                                             dataReplicates = 10,
                                             maxDispersion = 10^3))
#save(fit, file = "data analysis/results/boolean dispersed fit6.Robj")
load("data analysis/results/boolean dispersed fit3.Robj")
subsets <- unique(booldata$Subset)
subsetIndex <- 1:length(subsets)
#subsetIndex <- c(1, 21, 13, 19)
subsets <- unique(booldata$Subset)[subsetIndex]

require(pROC)
posteriors <- fit$posteriors[, -1, drop = FALSE]
par(mfrow = c(4, 6), mar = rep(1, 4))
#par(mfrow = c(2, 2), mar = rep(1, 4))
auc <- numeric(length(subsets))
for(j in 1:length(subsets)) {
  i <- which(names(posteriors) == subsets[j])
  try(rocfit <- roc(!vaccine ~ posteriors[, i]))
  auc[i] <- rocfit$auc
  print(plot(rocfit, main = paste(subsets[j], "- AUC", round(rocfit$auc, 3)),
             cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.6))
}

# Subject level posterior aggeregate?
par(mfrow = c(1, 1), mar = rep(5, 5))
weights <- apply(posteriors, 2, var)
weights <- weights / sum(weights)
aggregate <- as.vector(as.matrix(posteriors) %*% weights)
aggregate <- apply(posteriors, 1, function(x) sum(log(x + 10^-4)))
rocfit <- roc(vaccine ~ aggregate)
plot(pROC::roc(vaccine ~ aggregate),
     main = paste("AUC - Overall", round(rocfit$auc, 3)),
     cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.6)


par(mfrow = c(4, 6), mar = rep(2, 4))
#par(mfrow = c(2, 2), mar = rep(4, 4))
for(j in 1:length(subsets)) {
  i <- which(names(posteriors) == subsets[j])
  post <- posteriors[, i]
  treatment <- vaccine[order(post)]
  uniquePost <- sort(unique(post))
  nominalFDR <- sapply(uniquePost, function(x) mean(post[post <= x]))
  empFDR <- sapply(uniquePost, function(x) 1 - mean(vaccine[post <= x]))
  power <- sapply(uniquePost, function(x) sum(vaccine[post <= x]) / sum(vaccine))
  print(plot(nominalFDR, empFDR, type = "l", xlim = c(0, 1), ylim = c(0, 1), col = "red",
             main = paste(subsets[j])))
  lines(nominalFDR, power, col = "blue", lty = 2)
  abline(a = 0, b = 1)
  abline(v = c(0.05, 0.1), h = c(0.8, 0.9), col = "grey")
}

forplot <- list()
booldata <- booldata[order(as.character(booldata$PTID)), ]
for(j in 1:length(subsets)) {
  i <- which(names(posteriors) == subsets[j])
  post <- 1 - posteriors[, i]
  negprop <- log(booldata$count / booldata$parentcount)[booldata$Subset == subsets[j] & booldata$stim == "nonstim"]
  envprop <- log(booldata$count / booldata$parentcount)[booldata$Subset == subsets[j] & booldata$stim == "stim"]
  forplot[[j]] <- data.frame(subset = subsets[j],
                             negprop = negprop, envprop = envprop,
                             posterior = post, vaccine = vaccine)
}

#forplot <- do.call("rbind", forplot[c(1:21, 23)])
forplot <- do.call("rbind", forplot)
require(ggplot2)
print(ggplot(forplot) +
        geom_point(aes(x = negprop, y = envprop, col = posterior, shape = vaccine == 1)) +
        facet_wrap(~ subset, scales = 'free', ncol = 6) +
        geom_abline(slope = 1, intercept = 0) +
        theme_bw() + scale_colour_gradientn(colours=rainbow(4)))


library(gridExtra)
table <- data.frame(index = 1:length(subsets), subset = subsets,
                    responseProb = round(fit$levelProbs, 2),
                    AUC = round(auc, 2))
pdf("figures/cell subset table B.pdf", height = 8, width = 6)
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
nodes$label <- rownames(fit$isingfit$weiadj)
nodes$nfunctions <- sapply(regmatches(nodes$label, gregexpr(",", nodes$label)), length) + 1
names(edges)[5] <- "Dependence"
nodes$CD154 <- as.numeric(grepl("CD154", nodes$label))
ggplot() +
  scale_colour_gradient(limits=c(0, 1), low="white", high = "black") +
  geom_segment(data = edges, aes(x = xstart, y = ystart,
                                 xend = xend, yend = yend,
                                 alpha = Dependence,
                                 col = Dependence),
               size = 1) +
  scale_fill_gradient(low = "grey", high = "red", limits = c(0.48, 1)) +
  #scale_fill_gradient(low = "grey", high = "red", limits = c(1, 5)) +
  #scale_fill_gradient(low = "grey", high = "red", limits = c(0, 1)) +
  geom_point(data = nodes[1:5, ], aes(x = x, y = y, fill = auc), shape = 21,
             size = 12, col = "white") +
  geom_point(data = nodes[6:23, ], aes(x = x, y = y, fill = auc), shape = 21,
             size = 12, col = "white") +
  scale_size(range = c(0.3, 1)) +
  scale_shape(solid = FALSE) +
  geom_text(data = nodes, aes(x = x, y = y, label = nodes$label), size = 1.8) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(), legend.position = "none")

