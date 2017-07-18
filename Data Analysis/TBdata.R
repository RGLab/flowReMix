geomean <- function(x) exp(mean(log(x)))
library(flowReMix)
# Loading data -------------------------
tbdat <- readRDS("data/tb_rozot_booleans.RDS")
names(tbdat) <- tolower(names(tbdat))
tbdat$ptid <- sapply(strsplit(tbdat$samplename, "_"), function(x) x[[1]])
tbdat$count[is.na(tbdat$count)] <- 0

# Parsing Expressions ----------------
expressions <- unique(tbdat$pop)
newexp <- character(length(expressions))
for(i in 1:length(expressions)) {
  e <- expressions[i]
  e <- strsplit(e, "&")[[1]]
  res <- c()
  for(j in 1:length(e)) {
    sube <- strsplit(e[j], "-")[[1]][1]
    if(substr(sube, 1, 1) != "!"){
      res <- c(res, sube)
    }
  }
  newexp[i] <- paste(paste(res, collapse = "+"), "+", sep = "")
}
map <- cbind(expressions, newexp)
for(i in 1:nrow(map)) {
  tbdat$population[tbdat$population == map[i, 1]] <- map[i, 2]
}

# Defining subsets w stim ---------------------
# tbdat$subset <- paste(tbdat$parent, tbdat$population, sep = "/")
# tbdat <- subset(tbdat, stim != "EBV")
# tbdat$stimtemp <- tbdat$stim
# tbdat$stim[tbdat$stim %in% c("P1", "P2", "P3")] <- "P"
# datlist <- list()
# stims <- unique(tbdat$stim)
# stims <- stims[stims != "UNS"]
# for(i in 1:length(stims)) {
#   temp <- subset(tbdat, stim %in% c(stims[i], "UNS"))
#   temp$stim[temp$stim == "UNS"] <- "ctrl"
#   temp$stimgroup <- stims[i]
#   datlist[[i]] <- temp
# }
# tbdat <- do.call("rbind", datlist)
# tbdat$subset <- paste(tbdat$stimgroup, tbdat$subset, sep = "/")
# nonzerocounts <- by(tbdat, tbdat$subset, function(x) mean(x$count > 0))
# nonzerocounts <- data.frame(names(nonzerocounts), as.numeric(unlist(nonzerocounts)))

# Defining subsets wo stim --------------
tbdat$subset <- paste(tbdat$parent, tbdat$population, sep = "/")
tbdat <- subset(tbdat, stim != "EBV")
tbdat$stim[tbdat$stim == "UNS"] <- "ctrl"
tbdat$stim <- factor(tbdat$stim, levels = c("ctrl", "MP", "Mtbaux", "P1", "P2", "P3"))
nonzerocounts <- by(tbdat, tbdat$subset, function(x) mean(x$count > 0))
nonzerocounts <- data.frame(names(nonzerocounts), as.numeric(unlist(nonzerocounts)))

# Choosing subset of data for analysis -----------------
countkeep <- nonzerocounts[nonzerocounts[, 2] >= 0.2, 1]
popkeep <- c("cd4", "MAIT", "NKrainbow", "DCs", "Bcells", "conCD*")
stimkeep <- c("P", "MP", "Mtbaux")
# tempdat <- subset(tbdat, (subset %in% countkeep) & (parent %in% popkeep) & (stimgroup %in% stimkeep))
tempdat <- subset(tbdat, (subset %in% countkeep) & (parent %in% popkeep))
#tempdat$stim <- factor(tempdat$stim, levels = c("ctrl", stimkeep))
tempdat$subset <- factor(tempdat$subset)
rm(tbdat)

# Visuaizations ------------------------------
library(ggplot2)
tempdat$prop <- tempdat$count / tempdat$parentcount #+ 1 / tempdat$parentcount
subctrl <- subset(tempdat, stim == "ctrl")
substim <- subset(tempdat, stim != "ctrl")
wide <- merge(subctrl, substim, by = c("ptid", "subset"))
# pdf("figures/TBscat.pdf", width = 15, height = 15)
# ggplot(wide) +
#   geom_point(aes(x = prop.x, y = prop.y, col = type.x, shape = type.x)) +
#   facet_wrap(~ subset, scales = "free", nrow = floor(sqrt(length(unique(wide$subset))))) +
#   geom_abline(intercept = 0, slope = 1) +
#   theme_bw()
# dev.off()


# Analysis ----------------------------------
npost <- 8
control <- flowReMix_control(updateLag = 15, nsamp = 50, initMHcoef = 1,
                             nPosteriors = npost, centerCovariance = TRUE,
                             maxDispersion = 10^3, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             lastSample = round(40 / npost), isingInit = -log(99),
                             initMethod = "sparse")

# tempdat$stim <- tempdat$stimtemp
# tempdat$stim[tempdat$stim == "UNS"] <- "aUNS"
# tempdat$stim <- factor(tempdat$stim, levels = sort(unique(tempdat$stim)))
fit <- flowReMix(cbind(count, parentcount - count) ~ stim,
                 subject_id = ptid,
                 cell_type = subset,
                 cluster_variable = stim,
                 data = tempdat,
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "sparse",
                 iterations = 25,
                 parallel = TRUE,
                 verbose = TRUE, control = control)
#save(fit, file = "data analysis/results/TB all 6 nostim.Robj")
load(file = "data analysis/results/TB cluster 2.Robj")

# Scatter plots with posteriors ---------------
library(reshape2)
library(ggplot2)
post <- fit$posteriors
post <- melt(post, id = c("ptid"))
names(post) <- c("ptid", "subset", "posterior")
post <- merge(wide, post, by = c("ptid", "subset"))
# pdf("figures/TBscatPost.pdf", width = 15, height = 15)
# ggplot(post) +
#   geom_point(aes(x = prop.x, y = prop.y, col = 1 - posterior, shape = type.x)) +
#   facet_wrap(~ subset, scales = "free", nrow = floor(sqrt(length(unique(wide$subset))))) +
#   geom_abline(intercept = 0, slope = 1) +
#   theme_bw() +
#   scale_colour_gradientn(colours = rainbow(4))
# dev.off()

# Fitting ROC curves -----------------
post <- fit$posteriors
outcome <- by(tempdat, tempdat$ptid, function(x) unique(x$type))
outcome <- data.frame(ptid = names(outcome), outcome = as.character(outcome))
post <- merge(post, outcome)
outcome <- post[, ncol(post)]
post <- post[, -c(1, ncol(post))]
library(pROC)
aucs <- numeric(ncol(post))
for(i in 1:ncol(post)) {
  aucs[i] <- roc(outcome ~ post[, i])$auc
}
coefs <- sapply(fit$coefficients, function(x) x[2])

lprobs <- 1 - colMeans((fit$posteriors[, -1]))
keep <- lprobs > 0.2 & (apply(fit$posteriors[, -1], 2, min) > 0.1)
keep <- apply(fit$posteriors[, -1], 2, function(x) cummean(sort(x))[3]) < 0.1

n0 <- sum(outcome == "TB")
n1 <- sum(outcome == "LTBI")
pvals <- pwilcox(aucs * n0 * n1, n0, n1, lower.tail = FALSE)
sort(p.adjust(pvals[keep], method = "BH"))
rocResults <- data.frame(subsets = colnames(fit$posteriors[-1]),
                         auc = aucs, pvals = pvals, prob = fit$levelProbs)
rocResults$qvals <- NA
rocResults$qvals[keep] <- p.adjust(pvals[keep], method = "BH")
rocResults$coef <- sapply(fit$coefficients, function(x) max(x[-1]))
rocResults[order(rocResults$auc, decreasing = TRUE), ]
subtable <- subset(rocResults, !is.na(qvals))
subtable[order(subtable$auc, decreasing = TRUE), ]

# Permuting ROC --------------
keep <- rep(TRUE, nrow(rocResults))
reps <- 200
aucmat <- matrix(nrow = reps, ncol = sum(keep))
post <- fit$posteriors[, -1][keep]
n0 <- sum(outcome == "TB")
n1 <- sum(outcome == "LTBI")
perm <- c(rep(0, n0), rep(1, n1))
for(i in 1:reps) {
  perm <- perm[order(runif(n0 + n1))]
  aucmat[i, ] <- apply(post, 2, function(x) roc(perm ~ x)$auc)
  cat(i, " ")
}
cat("\n")
maxauc <- apply(aucmat, 1, mean)
hist(maxauc)
abline(v = mean(aucs))
mean(maxauc > mean(aucs))
quant <- quantile(maxauc, 0.9)

# Graph -----------------
assignments <- fit$assignmentList
names(assignments) <- sapply(names(assignments), function(x) strsplit(x, "%%%")[[1]][[1]])
assignments <- lapply(unique(names(assignments)), function(x) {
  do.call("rbind", assignments[names(assignments) == x])
})

subsets <- names(fit$coefficients)
reps <- 50
modelList <- list()
nsubsets <- ncol(assignments[[1]])
countCovar <- matrix(0, nrow = nsubsets, ncol = nsubsets)
doParallel::registerDoParallel(cores = 2)
for(i in 1:reps) {
  mat <- t(sapply(assignments, function(x) x[sample(1:nrow(x), 1), ]))
  colnames(mat) <- subsets
  keep <- apply(mat, 2, function(x) any(x != x[1]))
  mat <- mat[, keep]
  #system.time(model <- IsingFit::IsingFit(mat, AND = TRUE, plot = TRUE))
  system.time(coefs <- raIsing(mat, AND = TRUE, gamma = 0.25, method = "sparse",
                               cv = FALSE))
  #plot(model)
  #countCovar[keep, keep] <- countCovar[keep, keep] + (model$weiadj != 0) * sign(model$weiadj)
  countCovar[keep, keep] <- countCovar[keep, keep] + (coefs != 0) * sign(coefs)
  print(i)
}
doParallel::stopImplicitCluster()
props <- countCovar / reps
diag(props) <- 0
#save(props, file = "Data Analysis/results/TB graph 2.Robj")
table(props) / 2
threshold <- 0.5
which(abs(props) > threshold, arr.ind = TRUE)
props[abs(props) <= threshold] <- 0
sum(props != 0) / 2
#load(file = "Data Analysis/results/HVTN bool robust graph.Robj")

require(GGally)
library(network)
library(sna)
net <- props
# net <- raIsing(do.call("rbind", fit$assignmentList), AND = FALSE,
#                gamma = 0.25)
diag(net) <- 0
network <- net
colnames(net) <- subsets
rownames(net) <- subsets
#threshold <- 0.49
keep <- apply(network, 1, function(x) any(abs(x) > threshold))
network <- network[keep, keep]
net <- network(net)
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
aucs <- aucs
nodes$auc <- aucs[keep]

names(edges)[5] <- "Dependence"
lims <- max(abs(props))
library(ggplot2)
ggplot() +
  scale_colour_gradient2(low="dark red", high = "dark green") +
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

library(igraph)
network[network != 0] <- 1
graph <- graph.adjacency(network)
comp <- components(graph)

library(igraph)
network[network != 0] <- 1
graph <- graph.adjacency(network)
comp <- components(graph)

# ROCs for vaccination
n0 <- sum(outcome == "TB")
n1 <- sum(outcome == "LTBI")
par(mfrow = c(4, 4))
cpvals <- c()
groups <- list()
for(i in 1:sum(comp$csize >= 2)) {
  group <- which(comp$csize >= 2)[i]
  group <- which(comp$membership == group)
  group <- nodes$label[group]
  groups <- c(groups, list(group))
  cols <- which(names(fit$posteriors)[-1] %in% group)
  score <- apply(fit$posteriors[, -1][, cols], 1, geomean)
  rocfit <- roc(outcome ~ score)
  auc <- rocfit$auc
  pval <- pwilcox(auc * n0 * n1, n0, n1, lower.tail = FALSE)
  cpvals <- c(cpvals, pval)
  plot(rocfit, main = paste("Size", length(group), "-", round(auc, 3), "-", round(pval, 4)))
}

cqvals <- p.adjust(cpvals, method = "BH")
groups[cqvals < 0.1]


# Cumulative Reponse curves
assignments <- fit$assignmentList
names(assignments) <- sapply(names(assignments), function(x) strsplit(x, "%%%")[[1]][[1]])
assignments <- lapply(unique(names(assignments)), function(x) {
  do.call("rbind", assignments[names(assignments) == x])
})
post <- data.frame(t(sapply(assignments, colMeans)))
post <- cbind(1:nrow(post), post)
keep <- rep(rep(TRUE), ncol(assignments[[1]]))
assignments <- lapply(assignments, function(x) x[, keep])
subsets <- names(fit$coefficients)[keep]

nfunctions <- sapply(subsets, function(x) length(strsplit(x, "+", fixed = TRUE)[[1]]))
weightList <- list()
#weightList$poly <- weights <- nfunctions / choose(5, nfunctions)
weightList$func <- rep(1, length(nfunctions))

resultList <- list()
type <- "poly"
for(w in 1:length(weightList)) {
  weights <- weightList[[w]]
  weightname <- names(weightList)[w]
  resultList[[w]] <- list()
  print(weightname)
  for(i in 1:(length(assignments) + 2)) {
    cat(i, " ")
    if(i <= length(assignments)) {
      samp <- assignments[[i]]
      ptid <- fit$posteriors$ptid[i]
    } else if(i == (length(assignments) + 1)) {
      ptid <- "TB"
      samp <- do.call("rbind", assignments[outcome == "TB"])
    } else if(i == (length(assignments) + 2)){
      ptid <- "LTBI"
      samp <- do.call("rbind", assignments[outcome == "LTBI"])
    }
    colnames(samp) <- sapply(subsets, function(x) strsplit(x, "/")[[1]][[1]]) # for stim groups
    #colnames(samp) <- sapply(subsets, function(x) strsplit(x, "/")[[1]][[2]]) # for parent group
    #colnames(samp) <- rep("all", ncol(samp)) # for just one plot
    groups <- unique(colnames(samp))
    subjDatList <- list()
    for(j in 1:length(groups)) {
      group <- groups[j]
      subsamp <- samp[, colnames(samp) == group]
      subw <- weights[colnames(samp) == group]
      subw <- subw / sum(subw)
      groupSize <- ncol(subsamp)
      if(length(subw) == 1) {
        values <- subsamp * subw
      } else {
        values <- subsamp %*% subw
      }
      uniqueValues <- c(0, sort(unique(values)), 1)
      props <- sapply(uniqueValues, function(x) mean(values >= x))

      subjDatList[[j]] <- data.frame(ptid = ptid, group = group,
                                     index = weightname,
                                     presponses = uniqueValues,
                                     postProb = props)
    }

    resultList[[w]][[i]] <- do.call("rbind", subjDatList)
  }
  cat("\n")
}
responsedat <- do.call("rbind", do.call("rbind", resultList))
tboutcome <- by(tempdat, tempdat$ptid, function(x) x$type[1])
tboutcome <- data.frame(ptid = names(tboutcome), type = as.character(tboutcome))
forplot <- merge(responsedat, tboutcome, by.x = "ptid", by.y = "ptid",
                 all.x = TRUE)
summarized <- subset(forplot, ptid %in% c("TB", "LTBI"))
forplot <- subset(forplot, !(ptid %in%  c("TB", "LTBI")))
forplot <- forplot[order(forplot$ptid, forplot$group, forplot$presponses), ]
ggplot(forplot) + geom_step(aes(x = presponses, y = 1 - postProb, group = ptid,
                             col = factor(type)),
                         alpha = 0.42) +
  geom_step(data = summarized, aes(x = presponses, y = 1 - postProb,
                                linetype = ptid), size = 0.7) +
  facet_wrap( ~ group) + theme_bw() +
  xlab("At Least % Responsive Subsets") +
  ylab("Posterior Probabilities") #+

ksresult <- by(summarized, summarized$group, function(x) ks.test(x$presponse[x$ptid == "TB"], x$presponse[x$ptid == "LTBI"]))




