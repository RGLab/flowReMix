assign <- function(x) {
  x$prop <- x$count / x$parentcount
  assign <- as.numeric(by(x, x$subset, function(y) y$prop[1] > y$prop[2]))
  assign[assign == 1] <- -1
  result <- data.frame(ptid = x$ptid[1], subset = unique(x$subset), assign = assign)
  return(result)
}

require(pROC)
require(reshape2)
load("data/rv144_booleans.rda")
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
# booldata <- merge(booldata, data.frame(ptid = correlates$ptid,
#                                        IgAprim = correlates$IgAprim,
#                                        V2prim = correlates$V2prim))
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
names(correlates)[1] <- "ptid"
vaccine <- correlates[, c(1, 62, 64)]
vaccine$ptid <- factor(vaccine$ptid, levels = levels(booldata$ptid))
booldata <- with(booldata, booldata[order(subset, ptid, stim, decreasing = FALSE), ])
booldata <- merge(booldata, vaccine, sort = FALSE)
booldata$hiv <- NA
booldata$hiv[booldata$infect.y == "INFECTED"] <- TRUE
booldata$hiv[booldata$infect.y == "NON-INFECTED"] <- FALSE
booldata <- with(booldata, booldata[order(subset, ptid, stim, decreasing = FALSE), ])

# Analysis -------------
library(flowReMix)
npost <- 1
niter <- 30
control <- flowReMix_control(updateLag = 5, nsamp = 20, initMHcoef = 2.5,
                             keepEach = 5,
                             nPosteriors = npost, centerCovariance = TRUE,
                             maxDispersion = 1000, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             lastSample = 4, isingInit = -log(99),
                             ncores = 2,
                             preAssignCoefs = 1,
                             prior = 1, isingWprior = FALSE,
                             markovChainEM = FALSE,
                             initMethod = "robust",
                             learningRate = 0.6, keepWeightPercent = 0.9)

booldata <- subset(booldata, subset %in% c("TNFa,IFNg,IL4,IL2,CD154", "IFNg,IL4,IL2,CD154", "IL2,CD154", "TNFa,IL2,CD154"))
# booldata <- subset(booldata, subset %in% c("TNFa,IL2,CD154"))
booldata$subset <- factor(booldata$subset)
preAssignment <- do.call("rbind", by(booldata, booldata$ptid, assign))
booldata$stim <- factor(booldata$stim, levels = c("nonstim", "stim"))
system.time(fit <- flowReMix(cbind(count, parentcount - count) ~ stim,
                             subject_id = ptid,
                             cell_type = subset,
                             cluster_variable = stim,
                             data = booldata,
                             covariance = "sparse",
                             ising_model = "sparse",
                             regression_method = "robust",
                             iterations =  200,
                             cluster_assignment = preAssignment,
                             parallel = FALSE,
                             verbose = TRUE, control = control))
# save(fit, file = "data analysis/results/local_rv144_prior.Robj")
# plot(fit, type = "scatter")

fit$data <- booldata
# Plots -------------
scatter <- plot(fit, type = "scatter", target = vaccine)
rocplot <- plot(fit, type = "ROC", target = vaccine, direction = "<")
fdrplot <- plot(fit, type = "FDR", target = vaccine)
# rocplot
# scatter

# ROC for vaccinations -----------------------------
# sink("data analysis/results/RV144logisticSummary.txt")
ids <- fit$posteriors[, 1:2]
vaccine[, 1] <- as.character(vaccine[, 1])
vaccine[, 1] <- factor(vaccine[, 1], levels = levels(ids[, 1]))
vaccine <- vaccine[!is.na(vaccine[, 1]), ]
vaccine <- vaccine[order(vaccine[, 1]), ]
ids <- merge(ids, vaccine, all.x = TRUE, all.y = FALSE,
                 by = "ptid", sort = FALSE)
vaccination <- ids[, 3]

rocplot <- plot(fit, target = vaccine, type = "ROC", ncols = 6,
                direction = "auto", thresholdPalette = NULL,
                subsets = NULL)
# save_plot("figures/cowROCplot.pdf", rocplot,
#           base_height = 6,
#           base_width = 12)

rocResults <- summary(fit, target = vaccine, direction = ">", adjust = "BH",
                       sortAUC = FALSE)
rocResults[order(rocResults$auc, decreasing = TRUE), ]

# ROC for infection status -------------------
infectDat <- data.frame(ptid = rv144_correlates_data$PTID, infect = rv144_correlates_data$infect.y)
datId <- as.character(fit$posteriors$ptid)
infectID <- as.character(infectDat$ptid)
infectDat <- infectDat[infectID %in% datId, ]
infectDat$ptid <- factor(as.character(infectDat$ptid), levels = levels(booldata$ptid))
infectDat <- infectDat[order(infectDat$ptid), ]
ids <- merge(ids, infectDat, sort = FALSE)
infect <- ids[, 4]
infect[infect == "PLACEBO"] <- NA
infect <- factor(as.character(infect), levels = c("INFECTED", "NON-INFECTED"))

infectResults <- summary(fit, target = hiv, direction = "auto", adjust = "BH",
                          sortAUC = FALSE, pvalue = "wilcoxon")
infectResults$responseProb <- colMeans(fit$posteriors[, -1])
infectResults[order(infectResults$pvalue, decreasing = FALSE), ]

stab <- stabilityGraph(fit, type = "ising", cv = FALSE, reps = 100, cpus = 2,
                       gamma = 0.25, AND = TRUE)
# save(stab, file = "data analysis/results/rv144_15_niter30npost6_stab.Robj")
# Graph
threshold <- 0.85
load(file = "data analysis/results/rv144_15_niter30npost2_stab.Robj")
plot(stab, fill = rocResults$auc, threshold = threshold, seed = 1)
load(file = "data analysis/results/rv144_15_niter45npost2_stab.Robj")
plot(stab, fill = rocResults$auc, threshold = threshold, seed = 1)
load(file = "data analysis/results/rv144_15_niter60npost2_stab.Robj")
plot(stab, fill = rocResults$auc, threshold = threshold, seed = 1)
load(file = "data analysis/results/rv144_15_niter45npost4_stab.Robj")
plot(stab, fill = rocResults$auc, threshold = threshold, seed = 1)
load(file = "data analysis/results/rv144_15_niter60npost4_stab.Robj")
plot(stab, fill = rocResults$auc, threshold = threshold, seed = 1)
load(file = "data analysis/results/rv144_15_niter30npost6_stab.Robj")
plot(stab, fill = rocResults$auc, threshold = 0.9, seed = 1)
load(file = "data analysis/results/rv144_16_niter30npost4_stab.Robj")
plot(stab, fill = rocResults$auc, threshold = 0.85, seed = 1)
load(file = "data analysis/results/rv144_16_niter30npost6_stab.Robj")
plot(stab, fill = rocResults$auc, threshold = 0.82, seed = 1)
load(file = "data analysis/results/rv144_16_niter60npost4_stab.Robj")
plot(stab, fill = rocResults$auc, threshold = 0.85, seed = 1)

load(file = "data analysis/results/rv144_18_niter50npost2_stab.Robj")
plot(stab, fill = rocResults$auc, threshold = 0.85, seed = 1)
load(file = "data analysis/results/rv144_18_niter50npost4_stab.Robj")
plot(stab, fill = rocResults$auc, threshold = 0.85, seed = 1)
load(file = "data analysis/results/rv144_18_niter50npost6_stab.Robj")
plot(stab, fill = rocResults$auc, threshold = 0.85, seed = 1)
load(file = "data analysis/results/rv144_18_niter100npost2_stab.Robj")
plot(stab, fill = rocResults$auc, threshold = 0.92, seed = 1)
load(file = "data analysis/results/rv144_18_niter100npost4_stab.Robj")
plot(stab, fill = rocResults$auc, threshold = 0.92, seed = 1)
load(file = "data analysis/results/rv144_18_niter100npost6_stab.Robj")
plot(stab, fill = rocResults$auc, threshold = 0.85, seed = 1)
load(file = "data analysis/results/rv144_18_niter200npost2_stab.Robj")
plot(stab, fill = rocResults$auc, threshold = 0.85, seed = 1)
load(file = "data analysis/results/rv144_18_niter200npost4_stab.Robj")
plot(stab, fill = rocResults$auc, threshold = 0.85, seed = 1)



#######################
func <- rowSums(fit$posteriors[, -1])
funcAUC <- roc(infect ~ func)$auc
n0 <- sum(infect == "INFECTED", na.rm = TRUE)
n1 <- sum(infect == "NON-INFECTED", na.rm = TRUE)
pwilcox(funcAUC * n0 * n1, n0, n1, lower.tail = FALSE)

nfunctions <- sapply(subsets, function(x) length(gregexpr(",", paste(",", x))[[1]]))
M <- 6
weights <- nfunctions / (choose(M, nfunctions))
poly <- apply(fit$posteriors[, -1], 1, function(x) weighted.mean(x, weights))
# poly <- apply(fit$posteriors[, -1], 1, function(x) median(weights * x))
polyAUC <- roc(infect ~ poly)$auc
n0 <- sum(infect == "INFECTED", na.rm = TRUE)
n1 <- sum(infect == "NON-INFECTED", na.rm = TRUE)
pwilcox(polyAUC * n0 * n1, n0, n1, lower.tail = FALSE)

plot(fit, type = "boxplot", target = hiv, test = "wilcoxon",
     one_sided = TRUE, groups = "all", jitter = TRUE)

roc(vaccination ~ func)$auc
roc(vaccination ~ poly)$auc

group <- c(24, 21, 15, 8)
score <- rowMeans(fit$posteriors[, group])
rocfit <- roc(infect ~ score)
pwilcox(rocfit$auc * n0 * n1, n0, n1, lower.tail = FALSE)
ids$groupscore <- score
ids$poly <- poly
ids$func <- func

vaccines <- subset(correlates, infect.y != "PLACEBO")
vaccines$ptid <- as.character(vaccines$ptid)
ids$ptid <- as.character(ids$ptid)
vaccines <- merge(vaccines, ids, all.x = TRUE, all.y = TRUE,
                  by.x = "ptid", by.y = "ptid")

vaccines <- subset(vaccines, vaccines$infect != "PLACEBO")

plot(vaccines$PFS, vaccines$poly)
lines(lowess(vaccines$PFS, vaccines$poly), col = "red", lwd = 2)
cor(vaccines$PFS, vaccines$poly)
abline(v = c(0.05, .085))
abline(h = c(0.35))
target <- vaccines$poly > 0.35 & vaccines$PFS > 0.05 & vaccines$PFS < 0.085

summary(glm(infect ~ func + IgAprim + risk.medium + risk.high + sex,
            family = "binomial",
            data = vaccines))

summary(glm(infect ~ poly + IgAprim + risk.medium + risk.high + sex,
            family = "binomial",
            data = vaccines))
polylogpval <- summary(glm(infect ~ poly + IgAprim + risk.medium + risk.high + sex,
                           family = "binomial",
                           data = vaccines))$coefficients[2, 4] / 2

summary(glm(infect ~ groupscore + IgAprim + risk.medium + risk.high + sex,
            family = "binomial",
            data = vaccines))

rocfit <- roc(vaccines$infect ~ vaccines$FS)
pwilcox(rocfit$auc * n0 * n1, n0, n1, lower.tail = FALSE)
rocfit <- roc(vaccines$infect ~ vaccines$PFS)
pwilcox(rocfit$auc * n0 * n1, n0, n1, lower.tail = FALSE)

# Logistic regressions --------------
vaccines <- subset(correlates, infect.y != "PLACEBO")
resultList <- list()
adjRocList <- list()
for(i in 1:(ncol(fit$posteriors) - 1)) {
  vaccines$post <- NULL
  post <- fit$posteriors[!is.na(infect), c(1, i + 1)]
  names(post)[2] <- "post"
  vaccines <- merge(vaccines, post, by.x = "ptid", by.y = "ptid", all.x = TRUE)
  resultList[[i]] <- summary(glm(infect.y ~ post + IgAprim + agecat + risk.medium + risk.high + sex,
                                 family = "binomial",
                                 data = vaccines))
  resid <- lm(post ~ IgAprim + agecat + risk.medium + risk.high + sex,
              data = vaccines)$residuals
  infectResid <- glm(infect.y  ~ IgAprim + agecat + risk.medium + risk.high + sex,
                     family = "binomial", data = vaccines)$residuals

  adjRocList[[i]] <- roc(vaccines$infect.y ~ resid)
}

names(resultList) <- colnames(fit$posteriors)[-1]
names(adjRocList) <- colnames(fit$posteriors)[-1]
regResult <- t(sapply(resultList, function(x) x$coefficient[2, c(1,4)]))
regResult <- data.frame(regResult)
regResult$auc <- sapply(adjRocList, function(x) x$auc)
regResult$aucPval <- pwilcox(regResult$auc * n0 * n1, n0, n1, lower.tail = FALSE)
regResult$aucQval <- p.adjust(regResult$aucPval, method = "BH")
regResult[order(regResult[, 2], decreasing = FALSE), ]

# Stability Graph --------------------------
load(file = "data analysis/results/rv144AggreageStability2.Robj")
load(file = "data analysis/results/rv144AggStab2.Robj")
stab <- stability
stabPlot <- plot(stab, fill = rocResults$auc, threshold = .9)
stabPlot
plot(stab, fill = infectResults$auc, threshold = .95)
save_plot(stabPlot, filename = "figures/rv144IsingStb2hard.pdf",
          base_width = 7, base_height = 6)


# FDR Curves -------------
fdrTable <- fdrTable(fit, vaccination)
fdrplot <- plot(fit, target = vaccination, type = "FDR")
# save_plot("figures/RV144FDRplot.pdf", fdrplot,
#           base_height = 9,
#           base_width = 12)

# Scatter plots -----------------
scatter <- plot(fit, target = vaccination, type = "scatter")
# save_plot("figures/zeroRuleScatterRV144_2.pdf", scatter,
#           base_height = 9,
#           base_width = 12)

# Boxplots ------------------
# infect[is.na(infect)] <- "PLACEBO"
groups <- list()
groups[[1]] <- names(fit$posteriors)[-1]
names(groups) <- paste("p-value:", round(polylogpval, 4))
subsets <- colnames(fit$posteriors)[-1]
nfunctions <- sapply(subsets, function(x) length(strsplit(x, ",", fixed = TRUE)[[1]]))
weights <- list()
weights[[1]] <- nfunctions / choose(6, nfunctions)
names(weights) <- "Polyfunctionality"
box <- plot(fit, target = infect, type = "boxplot", groups = groups,
     weights = weights)
# save_plot("figures/RV144boxplotALLforBio.pdf", box,
#           base_width = 9, base_height = 5)


# Stability selection for graphical model ------------------------
# stability <- stabilityGraph(fit, type = "ising", cpus = 2, reps = 50,
#                             cv = FALSE, gamma = 0.25)
load("data analysis/results/RV144clusterIsing12.Robj")

isingplot <- plot(stability, fill = rocResults$auc, threshold = 0.78)
isingplot
# save_plot("figures/RV144isingplot.pdf", isingplot,
#           base_width = 9, base_height = 5)


# randStability <- stabilityGraph(fit, type = "randomEffects", cpus = 2, reps = 100,
#                                 cv = TRUE)
# load("data analysis/results/RV144clusterRandom6.Robj")
# randplot <- plot(randStability, fill = rocResults$auc, threshold = 0.96)
# randplot
# save_plot("figures/RV144randPlot.pdf", randplot,
#           base_width = 7, base_height = 5)

# Analysis with graph clusters --------------------
# groups <- getGraphComponents(stability, threshold = 0.96)
groups <- list(c("INFg,CD154", "CD154", "IFNG,IL2", "TNFa,CD154", "IL4,IL2,CD154",
                 "IFNg", "IL4,CD154"),
               c("IL2,CD154", "TNFa,IFNg,IL2,CD154",
                 "TNFa,IL4,IL2,CD154", "TNFa,IL2,CD154"))
groups$all <- c(groups[[1]], groups[[2]])
names(groups) <- c("Th2", "Th1", "all")
infect[infect == "PLACEBO"] <- NA
pvals <- numeric(2)
for(i in 1:length(groups)) {
  vaccines$groupscore <- NULL
  sub <- which(colnames(fit$posteriors) %in% groups[[i]])
  score <- apply(fit$posteriors[, sub], 1, function(x) weighted.mean(x, weights[[1]][sub - 1]))
  # score <- rowMeans(fit$posteriors[, sub])
  groupscore <- data.frame(PTID = fit$posteriors[, 1], groupscore = score)
  vaccines <- merge(vaccines, groupscore, all.x = TRUE)
  glmfit <- summary(glm(infect.y ~ groupscore + IgAprim + risk.medium + risk.high + sex,
                data = vaccines, family = "binomial"))
  pvals[i] <- glmfit$coefficients[2, 4] / 2
}
names(groups) <- paste(names(groups), "p-value:", round(pvals, 4))
graphbox <- plot(fit, type = "boxplot", target = infect, groups = groups,
     weights = weights, test = "none")
# save_plot(graphbox, file = "figures/RV144isingBoxplot.pdf",
#           base_height = 5, base_width = 9)

# Posterior probabilities for nresponses ----------------
graph <- fit$isingCov
thresholds <- diag(graph)
diag(graph) <- 0
assignment <- IsingSampler::IsingSampler(500, graph, thresholds)

assignments <- fit$assignmentList
subsets <- names(fit$coefficients)
resultList <- list()

subsets <- names(fit$posteriors)[-1]
selected <- sapply(fit$coefficients, function(x) x[2] > 0) & fit$levelProbs > 0.05
assignments <- lapply(assignments, function(x) x[, selected])
nfunctions <- sapply(subsets[selected], function(x) length(gregexpr(",", paste(",", x))[[1]]))
M <- 6
weights <- nfunctions / (choose(M, nfunctions))
#weights <- rep(1, length(selected))
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
  if(i == length(assignments) + 1) {
    next
    ptid <- "prior"
  } else {
    ptid <- names(assignments)[i]
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
outcome <- infectDat
forplot <- merge(responsedat, outcome, by.x = "ptid", by.y = "ptid",
                 all.x = TRUE)
forplot$VACCINE <- forplot$vaccine
# forplot <- merge(forplot, infectDat, by.x = "ptid", by.y = "ptid",
#                  all.x = TRUE)
forplot$INFECT <- forplot$infect

summarized <- summarize(group_by(forplot, group, presponses, INFECT),
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
#plot(rocfit, main = "ROC for Functionality Score - AUC = 0.9743")
subinfect <- as.numeric(by(forplot, forplot$ptid, function(x) x$infect[1]))
subinfect <- subinfect[subinfect != 3]
infectAUC <- roc(subinfect ~ intCDF[vaccine == 1])$auc
infectAUC
n1 <- sum(subinfect == 1)
n2 <- sum(subinfect == 2)
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
obsdiff <- min((curveA[, 4] - curveB[, 4]))
plot(density(diff))
abline(v = obsdiff)
mean(obsdiff > diff)

# Posteriors box plots ------------------------------------
require(dplyr)
require(reshape2)
outcome <- infectDat
posteriors <- fit$posteriors
selected <- sapply(fit$coefficients, function(x) x[2] > 0) & fit$levelProbs > 0.05
posteriors <- posteriors[, c(TRUE, selected)]
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

# Random Effect Covariance ---------------------
random <- fit$randomEffectSamp
reps <- 100
modelList <- list()
nsubsets <- ncol(assignments[[1]])
p <- ncol(random[[1]])
countCovar <- matrix(0, nrow = p , ncol = p )
doParallel::registerDoParallel(cores = 2)
for(i in 1:reps) {
  mat <- t(sapply(random, function(x) x[sample(1:nrow(x), 1), ]))
  colnames(mat) <- subsets
  model <- raIsing(mat, family = "gaussian", gamma = 0, cv = TRUE)
  modelList[[i]] <- model
  diag(model) <- 0
  countCovar <- countCovar + (model != 0) * sign(model)
  print(i)
}
doParallel::stopImplicitCluster()
props <- countCovar / reps
#save(props, file = "data analysis/results/RV144 RE net 1.Robj")
table(props)
estnet <- props
diag(estnet) <- 0
threshold <- 0.6
estnet[abs(estnet) < threshold] <- 0

network <- estnet
require(GGally)
library(network)
library(sna)
keep <- apply(network, 1, function(x) any(abs(x) >= threshold)) #| rocResults$qvals < 0.05
network <- network[keep, keep]
net <- network(network)
subsets <- colnames(fit$posteriors)[-1]
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
  scale_colour_gradient2(limits = c(-lims, lims), low="dark red", high = "dark green") +
  geom_segment(data = edges, aes(x = xstart, y = ystart,
                                 xend = xend, yend = yend,
                                 col = (Dependence),
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




