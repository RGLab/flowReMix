assign <- function(x) {
  x$prop <- x$count / x$parentcount
  assign <- as.numeric(by(x, x$subset, function(y) y$prop[1] > y$prop[2]))
  assign[assign == 1] <- -1
  result <- data.frame(ptid = x$ptid[1], subset = unique(x$subset), assign = assign)
  return(result)
}

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
names(booldata) <- tolower(names(booldata))

# Getting vaccine information --------------------
data("rv144")
rv144 <- rv144[order(rv144$ptid), ]
vaccine <- (by(rv144, rv144$ptid, function(x) x$vaccine[1] == "VACCINE"))
vaccine <- data.frame(ptid = names(vaccine), vaccine = as.numeric(vaccine))
vaccinemat <- vaccine[vaccine$ptid %in% booldata$ptid, ]

# Getting infection status
data("rv144_correlates_data")
correlates <- rv144_correlates_data
correlates <- correlates[order(as.character(correlates$PTID)), ]
infection <- correlates$infect.y

# Analysis -------------
library(flowReMix)
control <- flowReMix_control(updateLag = 6, nsamp = 100, initMHcoef = 2.5,
                             nPosteriors = 1, centerCovariance = TRUE,
                             maxDispersion = 10^3, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             initMethod = "robust")

booldata$subset <- factor(booldata$subset)
preAssignment <- do.call("rbind", by(booldata, booldata$ptid, assign))
system.time(fit <- flowReMix(cbind(count, parentcount - count) ~ treatment,
                 subject_id = ptid,
                 cell_type = subset,
                 cluster_variable = treatment,
                 data = booldata,
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "robust",
                 iterations = 12,
                 cluster_assignment = preAssignment,
                 parallel = TRUE,
                 verbose = TRUE, control = control))
#save(fit, file = "data analysis/results/boolean robust2 wPre.Robj")
load(file = "data analysis/results/boolean robust2.Robj")
# load(file = "data analysis/results/boolean robust2 wPre.Robj")
# load(file = "data analysis/results/boolean upfit4 w pre.Robj")
# load(file = "data analysis/results/boolean upfit3.Robj")



# ROC for vaccinations -----------------------------
subsets <- unique(booldata$subset)
subsetIndex <- 1:length(subsets)
#subsetIndex <- c(1, 21, 13, 19)
subsets <- unique(booldata$subset)[subsetIndex]

require(pROC)
posteriors <- fit$posteriors
posteriors <- posteriors[order(fit$posteriors$ptid), ]
vaccine <- vaccinemat
vaccine[, 1] <- factor(as.character(vaccine[, 1]), levels = levels(fit$posteriors$ptid))
vaccine <- vaccine[order(vaccine[, 1]), ]
vaccine <- vaccine[, 2]
par(mfrow = c(4, 6), mar = rep(1, 4))
#par(mfrow = c(2, 2), mar = rep(1, 4))
auc <- numeric(length(subsets))
for(j in 1:length(subsets)) {
  i <- which(names(posteriors) == subsets[j])
  try(rocfit <- roc(!vaccine ~ posteriors[, i]))
  auc[j] <- rocfit$auc
  # try(print(plot(rocfit, main = paste(subsets[j], "- AUC", round(rocfit$auc, 3)),
  #            cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.6)))
}

n1 <- sum(vaccine)
n2 <- sum(!vaccine)
wilcox <- auc * n1 * n2
pvals <- pwilcox(wilcox, n1, n2, lower.tail = FALSE)
rocResults <- data.frame(subsets, auc, pvals)
qvals <- p.adjust(pvals, method = "bonferroni")
rocResults$qvals <- qvals
rocResults[order(rocResults$pvals), ]
select <- rocResults$qvals < 0.05

# ROC analysis for infection status --------------------
require(pROC)
posteriors <- fit$posteriors[, -1, drop = FALSE]
posteriors <- subset(posteriors, infection != "PLACEBO")
subinfect <- infection[infection != "PLACEBO"]
par(mfrow = c(4, 6), mar = rep(1, 4))
#par(mfrow = c(2, 2), mar = rep(1, 4))
auc <- numeric(length(subsets))
for(j in 1:length(subsets)) {
  i <- which(names(posteriors) == subsets[j])
  try(rocfit <- roc(subinfect ~ posteriors[, i]))
  auc[i] <- rocfit$auc
  # try(print(plot(rocfit, main = paste(subsets[j], "- AUC", round(rocfit$auc, 3)),
  #                cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.6)))
}

n1 <- sum(subinfect == "INFECTED")
n2 <- length(subinfect) - n1
wilcox <- auc * n1 * n2
pvals <- pwilcox(wilcox, n1, n2, lower.tail = FALSE)
infectResult <- data.frame(subsets, auc, pvals)
infectResult <- subset(infectResult, select)
qvals <- p.adjust(infectResult$pvals, method = "BY")
infectResult$qvals <- qvals
infectResult[order(infectResult$pvals), ]


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


posteriors <- fit$posteriors
posteriors <- posteriors[order(fit$posteriors$ptid), ]
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

# Scatter plots -----------------
forplot <- list()
booldata <- booldata[order(as.character(booldata$ptid)), ]
fit$posteriors <- fit$posteriors[order(as.character(fit$posteriors$ptid)), ]
posteriors <- fit$posteriors[, -1, drop = FALSE]
logit <- function(x) log(x / (1 - x))
for(j in 1:length(subsets)) {
  i <- which(names(posteriors) == subsets[j])
  post <- 1 - posteriors[, i]
  negprop <- logit(booldata$count / booldata$parentcount)[booldata$subset == subsets[j] & booldata$stim == "nonstim"]
  envprop <- logit(booldata$count / booldata$parentcount)[booldata$subset == subsets[j] & booldata$stim == "stim"]
  forplot[[j]] <- data.frame(subset = subsets[j],
                             negprop = negprop, envprop = envprop,
                             posterior = post, vaccine = vaccine,
                             ptid = fit$posteriors$ptid)
}

#forplot <- do.call("rbind", forplot[c(1:21, 23)])
infectDat <- data.frame(ptid = rv144_correlates_data$PTID, infect = rv144_correlates_data$infect.y)
forplot <- do.call("rbind", forplot)
forplot <- merge(forplot, infectDat, all.x = TRUE, by.x = "ptid", by.y = "ptid")
require(ggplot2)
ggplot(forplot) +
        geom_point(aes(x = negprop, y = envprop, col = posterior, shape = vaccine == 1)) +
        facet_wrap(~ subset, scales = 'free', ncol = 6) +
        geom_abline(slope = 1, intercept = 0) +
        theme_bw() + scale_colour_gradientn(colours=rainbow(4))

# library(gridExtra)
# table <- data.frame(index = 1:length(subsets), subset = subsets,
#                     responseProb = round(fit$levelProbs, 2),
#                     AUC = round(auc, 2))
# pdf("figures/cell subset table B.pdf", height = 8, width = 6)
# grid.table(table, rows = NULL)
# dev.off()

# Testing Screening Procedure -------------------
assignments <- fit$assignmentList
names(assignments) <- substr(names(assignments), 1, 5)
assignments <- lapply(unique(names(assignments)), function(x) {
  do.call("rbind", assignments[names(assignments) == x])
})
subsets <- names(fit$coefficients)

nSubsets <- length(fit$coefficients)
nSubjects <- length(fit$assignmentList)
screenPvals <- numeric(nSubsets)
M <- numeric(nSubsets)
for(j in 1:nSubsets) {
  counts <- sapply(assignments, function(x) sum(x[, j]))
  N <- sapply(assignments, function(x) length(x[, j]))
  counts <- pmax(pmin(N - 1, counts), 1)
  stest <- mcScreenTest(counts, N, 2000)
  screenPvals[j] <- stest$pval
  print(c(j, stest$M, screenPvals[j]))
}
pvals <- pchisq(likratios, nSubjects - 1, lower.tail = FALSE)
cbind(rocResults, screenPvals)
which(p.adjust(pvals, method = "bonferroni")  > 0.01)


# Stability selection for graphical model ------------------------
assignments <- fit$assignmentList
names(assignments) <- substr(names(assignments), 1, 5)
assignments <- lapply(unique(names(assignments)), function(x) {
  do.call("rbind", assignments[names(assignments) == x])
})
subsets <- names(fit$coefficients)

reps <- 40
modelList <- list()
nsubsets <- ncol(assignments[[1]])
countCovar <- matrix(0, nrow = nsubsets, ncol = nsubsets)
for(i in 1:reps) {
  mat <- t(sapply(assignments, function(x) x[sample(1:nrow(x), 1), ]))
  colnames(mat) <- subsets
  keep <- apply(mat, 2, function(x) any(x != x[1]))
  mat <- mat[, keep]
  model <- IsingFit::IsingFit(mat, AND = TRUE, plot = FALSE)
  modelList[[i]] <- model
  #plot(model)
  countCovar[keep, keep] <- countCovar[keep, keep] + (model$weiadj != 0) * sign(model$weiadj)
  print(i)
}

props <- countCovar / reps
table(props)
threshold <- 0.5
which(props > threshold, arr.ind = TRUE)
props[abs(props) < threshold] <- 0
sum(props != 0) / 2
#save(props, file = "Data Analysis/results/RV144 1 graph.Robj")

# Plotting graph ---------------------
require(GGally)
library(network)
library(sna)
network <- props
keep <- apply(network, 1, function(x) any(abs(x) >= threshold)) #| rocResults$qvals < 0.05
#keep <-  rep(TRUE, nrow(props))
network <- network[keep, keep]
net <- network(props)
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
nodes$auc <- rocResults$auc[keep]
nodes$qvals <- rocResults$qvals[keep]
nodes$sig <- nodes$qvals < 0.05

names(edges)[5] <- "Dependence"
lims <- max(abs(props))
library(ggplot2)
ggplot() +
  scale_colour_gradient2(limits=c(-lims, lims), low="dark red", high = "dark green") +
  geom_segment(data = edges, aes(x = xstart, y = ystart,
                                 xend = xend, yend = yend,
                                 col = abs(Dependence),
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

# Posterior probabilities for nresponses ----------------
graph <- fit$isingCov
thresholds <- diag(graph)
diag(graph) <- 0
assignment <- IsingSampler::IsingSampler(500, graph, thresholds)

assignments <- fit$assignmentList
names(assignments) <- substr(names(assignments), 1, 5)
assignments <- lapply(unique(names(assignments)), function(x) {
  do.call("rbind", assignments[names(assignments) == x])
})
selected <- rocResults$qvals < 0.05
assignments <- lapply(assignments, function(x) x[, selected])
subsets <- names(fit$coefficients)[selected]
resultList <- list()

nfunctions <- sapply(subsets, function(x) length(gregexpr(",", paste(",", x))[[1]]))
M <- 6
weights <- nfunctions / (choose(M, nfunctions))
weights <- rep(1, length(subsets))
weights <- weights / sum(weights)

allassign <- do.call("rbind", assignments)
values <- unique(as.numeric(allassign %*% weights))
values <- sort(c(values, 0, 1))
for(i in 1:length(assignments)) {
  cat(i, " ")
  samp <- assignments[[i]]
  colnames(samp) <- rep("all", ncol(samp)) # for just one plot
  groups <- unique(colnames(samp))
  subjDatList <- list()
  if(i == length(assignments)) {
    ptid <- "prior"
  } else {
    ptid <- fit$posteriors$ptid[i]
  }
  for(j in 1:length(groups)) {
    group <- groups[j]
    subsamp <- samp[, groups == group]
    w <- weights[groups == group]
    w <- w / sum(w)
    groupSize <- ncol(subsamp)
    scores <- as.numeric(subsamp %*% w)
    #values <- sort(c(0, unique(scores), 1))
    props <- sapply(values, function(x) mean(scores >= x))
    subjDatList[[j]] <- data.frame(ptid = ptid, group = group,
                                   presponses = values,
                                   postProb = props)
  }
  resultList[[i]] <- do.call("rbind", subjDatList)
}
responsedat <- do.call("rbind", resultList)
outcome <- data.frame(ptid = fit$posteriors$ptid, vaccine = vaccine)
forplot <- merge(responsedat, outcome, by.x = "ptid", by.y = "ptid",
                 all.x = TRUE)
forplot$VACCINE <- forplot$vaccine
forplot <- merge(forplot, infectDat, by.x = "ptid", by.y = "ptid",
                 all.x = TRUE)
forplot$INFECT <- forplot$infect

summarized <- summarize(group_by(forplot, group, presponses, VACCINE,
                                 INFECT),
                        postProb = mean(postProb))

ggplot(forplot) + geom_line(aes(x = presponses, y = 1 - postProb, group = ptid,
                                col = factor(infect)),
                            alpha = 0.42) +
  geom_line(data = summarized, aes(x = presponses, y = 1 - postProb,
                                   linetype = factor(INFECT)), size = 1) +
  facet_wrap(~ group) +
  xlab("At Least % Responsive Subsets") +
  ylab("Posterior Probabilities") + theme_bw()

integrateCDF <- function(x) {
  x <- x[, 3:4]
  result <- 0
  for(i in 2:nrow(x)) {
    base <- (x[i, 1] - x[i - 1, 1])
    result <- result + base * (x[i, 2])
    #result <- result + base * (x[i - 1, 2] - x[i, 2]) / 2
  }

  return(result)
}
intCDF <- by(forplot, forplot$ptid, integrateCDF)
intCDF <- as.numeric(intCDF)
roc(vaccine ~ intCDF)$auc
rocfit <- roc(vaccine ~ intCDF)
par(mfrow = c(1, 1))
plot(rocfit, main = "ROC for Functionality Score - AUC = 0.9743")
infectAUC <- roc(subinfect ~ intCDF[vaccine == 1])$auc
infectAUC
n1 <- sum(subinfect == "INFECTED")
n2 <- sum(subinfect == "NON-INFECTED")
pval <- pwilcox(infectAUC * n1 * n2, n1, n2, FALSE)
plot(roc(subinfect ~ intCDF[vaccine == 1]))

correlates$intCDF <- intCDF
vaccines <- subset(correlates, infect.y != "PLACEBO")
summary(glm(infect.y ~ intCDF + IgAprim + agecat + risk.medium + risk.high + sex,
    family = "binomial",
    data = vaccines))

lmfit <- lm(intCDF ~ IgAprim + agecat + risk.medium + risk.high + sex, data = vaccines)
plot(roc(vaccines$infect.y ~ lmfit$residuals)$auc)

# Bootstrapping `KS` test ---------------------
tempdat <- subset(forplot, infect != "PLACEBO")
reps <- 200
diff <- numeric(reps)
ids <- unique(forplot$ptid)
nsubjects <- length(ids)
infect <- infectDat[infectDat[, 2] != "PLACEBO", 2]
for(i in 1:reps) {
  tempinfect <- infect[order(runif(nsubjects))]
  for(j in 1:nsubjects) {
    tempdat$infect[tempdat$ptid == ids[j]] <- tempinfect[j]
  }

  curveA <- data.frame(summarize(group_by(subset(tempdat, infect == "INFECTED"), group, presponses, infect),
                      postProb = mean(postProb)))
  curveB <- data.frame(summarize(group_by(subset(tempdat, infect != "INFECTED"), group, presponses, infect),
                      postProb = mean(postProb)))
  diff[i] <- min(as.numeric((curveA[, 4] - curveB[, 4])))
  cat(i, " ")
}

curveA <- data.frame(subset(summarized, INFECT == "INFECTED"))
curveB <- data.frame(subset(summarized, INFECT == "NON-INFECTED"))
obsdiff <- min((curveA[, 5] - curveB[, 5]))
plot(density(diff))
abline(v = obsdiff)
mean(obsdiff > diff)

# Posteriors box plots ------------------------------------
require(dplyr)
require(reshape2)
outcome <- infectDat
posteriors <- fit$posteriors
posteriors[, c(16, 22)] <- 1 - posteriors[, c(16, 22)]
posteriors <- merge(posteriors, infectDat,
                    by.x = "ptid", by.y = "ptid",
                    all.x = TRUE)
posteriors <- melt(posteriors, id.vars = c("ptid", "infect"))
names(posteriors)[3:4] <- c("subset", "posterior")

ggplot(posteriors, aes(x = infect, y = 1- posterior)) +
  geom_boxplot() + geom_jitter(size = 0.05, alpha = 0.2) +
  xlab("vaccine") + theme_bw()

ggplot(posteriors, aes(x = infect, y = 1- posterior, col = infect)) +
  geom_boxplot() + #geom_jitter(size = 0.05, alpha = 0.2) +
  xlab("vaccine") + theme_bw() + facet_wrap(~ subset, nrow = 3)




