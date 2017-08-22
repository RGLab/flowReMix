getExpression <- function(str) {
  first <- substr(str, 1, 7)
  second <- substr(str, 8, nchar(str))
  second <- strsplit(second, "")[[1]]
  seperators <- c(0, which(second %in% c("-", "+")))
  expressed <- list()
  for(i in 2:length(seperators)) {
    if(second[seperators[i]] == "+") {
      expressed[[i]] <- paste(second[(seperators[(i - 1)] + 1) : seperators[i]], collapse = '')
    }
  }

  expressed <- paste(unlist(expressed), collapse = '')
  expressed <- paste(first, expressed, sep = '')
  return(expressed)
}

assign <- function(x) {
  x$prop <- x$count / x$parentcount
  assign <- as.numeric(by(x, x$subset, function(y) max(y$prop[y$stim != 0]) > min(y$prop[y$stim == 0])))
  assign[assign == 1] <- -1
  result <- data.frame(ptid = x$ptid[1], subset = unique(x$subset), assign = assign)
  return(result)
}

# Loading Data --------------------------------
# hvtn <- read.csv(file = "data/merged_505_stats.csv")
# names(hvtn) <- tolower(names(hvtn))
# hvtn <- subset(hvtn, !is.na(ptid))
# saveRDS(hvtn, file = "data/505_stats.rds")


# Getting Demographic data ------------------------
demo <- read.csv(file = "data/primary505.csv")
infect <- data.frame(ptid = demo$ptid, status = demo$HIVwk28preunbl)
infect <- subset(infect, infect$ptid %in% hvtn$ptid)

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
marginals <- subset(marginals, stim %in% c("negctrl", "VRC ENV A",
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
env <- subset(marginals, stim %in% c("negctrl", "VRC ENV C", "VRC ENV B", "VRC ENV A"))
env$subset <- factor(paste("env", env$population, sep = "/"))
env$stimGroup <- "env"
nef <- subset(marginals, stim %in% c("negctrl", "VRC NEF B"))
nef$subset <- factor(paste("nef", nef$population, sep = "/"))
nef$stimGroup <- "nef"
subsetDat <- rbind(gag, pol, env, nef)
subsetDat$stim <- as.character(subsetDat$stim)
subsetDat$stim[subsetDat$stim == "negctrl"] <- 0
subsetDat$stim <- factor(subsetDat$stim)

# Converting subset names ------------------
subsets <- as.character(unique(subsetDat$subset))
expressed <- sapply(subsets, getExpression)
map <- cbind(subsets, expressed)
subsetDat$subset <- as.character(subsetDat$subset)
for(i in 1:nrow(map)) {
  subsetDat$subset[which(subsetDat$subset == map[i, 1])] <- map[i, 2]
}
subsetDat$subset <- factor(subsetDat$subset)

# Getting outcomes -------------------------------
treatmentdat <- read.csv(file = "data/rx_v2.csv")
names(treatmentdat) <- tolower(names(treatmentdat))
treatmentdat$ptid <- factor(gsub("-", "", (treatmentdat$ptid)))
treatmentdat <- subset(treatmentdat, ptid %in% unique(subsetDat$ptid))

# Finding problematic subsets?
keep <- by(subsetDat, list(subsetDat$subset), function(x) mean(x$count > 1) > 0.02)
keep <- names(keep[sapply(keep, function(x) x)])
#result$subsets[result$qvals < 0.1] %in% keep
subsetDat <- subset(subsetDat, subset %in% keep)
subsetDat$subset <- factor(as.character(subsetDat$subset))

# Fitting the model ------------------------------
library(flowReMix)
control <- flowReMix_control(updateLag = 12, nsamp = 100, initMHcoef = 2.5,
                             nPosteriors = 1, centerCovariance = TRUE,
                             maxDispersion = 10^3, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             lastSample = 100, isingInit = -log(95),
                             initMethod = "robust")

subsetDat$batch <- factor(subsetDat$batch..)
subsetDat$stimGroup <- factor(subsetDat$stimGroup)
preAssign <- by(subsetDat, subsetDat$ptid, assign)
preAssign <- do.call("rbind", preAssign)
fit <- flowReMix(cbind(count, parentcount - count) ~ stim,
                 subject_id = ptid,
                 cell_type = subset,
                 cluster_variable = stim,
                 data = subsetDat,
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "robust",
                 iterations = 20,
                 parallel = FALSE,
                 verbose = TRUE, control = control)
# load(file = "Data Analysis/results/HVTN bool robust.Robj")
#load(file = "Data Analysis/results/HVTNclust1.Robj")
# load(file = "Data Analysis/results/HVTNclust2.Robj")
# load(file = "Data Analysis/results/HVTNclust4.Robj")
load(file = "Data Analysis/results/HVTNclust7npost1niter24.Robj")
load(file = "Data Analysis/results/HVTNclust7npost10niter24.Robj")
# fit$posteriors[, -1] <- 1 - fit$posteriors[, -1]

# Post-hoc assignment -------------------
# subjects <- unique(preAssign$ptid)
# for(i in 1:length(subjects)) {
#   row <- which(fit$posteriors$ptid == subjects[i])
#   assign <- subset(preAssign, ptid == subjects[i])
#   matching <- match(colnames(fit$posteriors[, -1]), assign[, 2])
#   index <- which(assign[matching, 3] == 0) + 1
#   fit$posteriors[row, index] <- fit$posteriors[row, index] / 100
# }


# ROC plots -----------------------------
require(pROC)
outcome <- treatmentdat[, c(13, 15)]
rocResults <- rocTable(fit, outcome[, 2], direction = ">", adjust = "BH",
                       sortAUC = FALSE)
rocResults[order(rocResults$auc, decreasing = TRUE), ]

vaccine <- outcome[, 2] == 0
hiv <- infect[, 2]
hiv[vaccine == 0] <- NA
infectROC <- rocTable(fit, hiv, direction = ">", adjust = "BH",
                      sortAUC = FALSE)
infectROC[order(infectROC$auc, decreasing = TRUE), ]

# Raw Graphical Models ---------------
isingThreshold <- 0.945
plot(fit, type = "graph", graph = "ising",
     fill = rocResults$auc, normalize = FALSE,
     threshold = isingThreshold)

plot(fit, type = "graph", graph = "ising",
     fill = infectROC$auc, normalize = FALSE,
     threshold = isingThreshold)

plot(fit, type = "graph", graph = "randomEffects",
     fill = rocResults$auc, normalize = FALSE,
     threshold = 0.99)


# Scatter plots -----------------------
scatter <- plot(fit, target = vaccine, type = "scatter", ncol = 11)
# save_plot("figures/HVTNscatter.pdf", scatter,
#           base_height = 20,
#           base_width = 22, limitsize = FALSE)

# FDR curves ----------------
vaccination <- outcome[, 2] == 0
fdrplot <- plot(fit, target = vaccination, type = "FDR")
# save_plot("figures/hvtnFDRplot.pdf", fdrplot,
#           base_height = 20,
#           base_width = 22, limitsize = FALSE)

# Posterior boxplots ------------------------
nfunctions <- sapply(strsplit(colnames(fit$posteriors)[-1], "+", fixed = TRUE), function(x) length(x) - 1)
weightList <- list()
weightList$Polyfunctionality <- weights <- nfunctions / choose(5, nfunctions)
# weightList$functionality <- rep(1, length(nfunctions))

subsets <- names(fit$posteriors[, -1])
stim <- sapply(strsplit(subsets, "/"), function(x) x[1])
stimnames <- unique(stim)
stim <- lapply(stimnames, function(x) subsets[stim %in% x])
names(stim) <- stimnames
stimpvals <- numeric(length(stim))
for(i in 1:length(stim)) {
  group <- stim[[i]]
  sub <- which(names(fit$posteriors)[-1] %in% group)
  aggregate <- apply(fit$posteriors[, -1], 1, function(x) weighted.mean(x[sub], weightList[[1]][sub]))
  stimpvals[i] <- summary(glm(hiv ~ aggregate, family = "binomial"))$coefficients[2, 4] / 2
}
names(stim) <- paste(names(stim), "pvalue:", round(stimpvals, 4))


parent <- sapply(strsplit(subsets, "/"), function(x) x[2])
parentnames <- unique(parent)
parent <- lapply(parentnames, function(x) subsets[parent %in% x])
names(parent) <- parentnames
parentpvals <- numeric(length(parent))
for(i in 1:length(parent)) {
  group <- parent[[i]]
  sub <- which(names(fit$posteriors)[-1] %in% group)
  aggregate <- apply(fit$posteriors[, -1], 1, function(x) weighted.mean(x[sub], weightList[[1]][sub]))
  parentpvals[i] <- summary(glm(hiv ~ aggregate, family = "binomial"))$coefficients[2, 4] / 2
}
names(parent) <- paste(names(parent), "pvalue:", round(parentpvals, 4))

stimparent  <- sapply(strsplit(subsets, "/"), function(x) paste(x[1:2], collapse = "/"))
stimparentnames <-unique(stimparent)
stimparent <- lapply(stimparentnames, function(x) subsets[stimparent %in% x])
names(stimparent) <- stimparentnames
sppvals  <- numeric(length(stimparent))
for(i in 1:length(stimparent)) {
  group <- stimparent[[i]]
  sub <- which(names(fit$posteriors)[-1] %in% group)
  aggregate <- apply(fit$posteriors[, -1], 1, function(x) weighted.mean(x[sub], weightList[[1]][sub]))
  sppvals[i] <- summary(glm(hiv ~ aggregate, family = "binomial"))$coefficients[2, 4] / 2
}
names(stimparent) <- paste(names(stimparent), "pvalue:", round(sppvals, 4))

infection <- hiv
infection[hiv == 0] <- "NON-INFECTED"
infection[hiv == 1] <- "INFECTED"
# infection[is.na(hiv)] <- "PLACEBO"
infection <- factor(infection, levels = c("PLACEBO", "INFECTED", "NON-INFECTED"))
stimbox <- plot(fit, type = "boxplot", groups = stim,
                weights = weightList, ncol = 2,
                target = infection)
stimbox
# save_plot("figures/hvtnStimBox.pdf", stimbox,
#           base_width = 9, base_height = 6)
parentbox <- plot(fit, type = "boxplot", groups = parent,
                  weights = weightList, ncol = 2,
                  target = infection)
parentbox
# save_plot("figures/hvtnParentBox.pdf", parentbox,
#           base_width = 9, base_height = 5)
parentstimbox <- plot(fit, type = "boxplot", groups = stimparent,
                      weights = weightList, ncol = 3, target = infection)
parentstimbox
# save_plot("figures/hvtnParentStimBox.pdf", parentstimbox,
#           base_width = 10, base_height = 8)

all <- list()
all[[1]] <- names(fit$posteriors)[-1]
poly <- apply(fit$posteriors[, -1], 1, function(x) weighted.mean(x, weightList[[1]]))
allpval <- summary(glm(hiv ~ poly, family = "binomial"))$coefficients[2, 4] / 2
names(all) <- paste("All Subsets, p-value:", round(allpval, 4))
allboxplot <- plot(fit, target = infection, type = "boxplot", groups = all)
# save_plot("figures/HVTNallboxplot.pdf", allboxplot)

# Stability selection for graphical model ------------------------
# load("data analysis/results/HVTNising2.Robj")
# stability <- stabilityGraph(fit, type = "ising", reps = 50, cpus = 2, gamma = 0.25)
# save(stability, file = "data analysis/results/HVTN bool robust graph4.Robj")
load("data analysis/results/HVTNising7setting3.Robj")
colnames(stability$network) <- colnames(fit$posteriors)[-1]
isingplot <- plot(stability,
                  fill = rocResults$auc,
                  threshold = 0.82)
isingplot
# save_plot("figures/HVTNising2.pdf", isingplot,
#           base_width = 9, base_height = 5)

load("data analysis/results/HVTNrand2.Robj")
randplot <- plot(randStability, fill = rocResults$auc,
                  threshold = 0.9)
randplot
# save_plot("figures/HVTNrand2.pdf", isingplot,
#           base_width = 7, base_height = 5)

# Boxplots for graph clusters -----------
groups <- getGraphComponents(stability, threshold = 0.65, minsize = 4)
weightList$Polyfunctionality <- weights <- nfunctions / choose(5, nfunctions)
names(groups) <- c("107+", "Large CD4+", "env/8+")
# pvals <- numeric(length(groups))
# weights <- list()
# weights[[1]] <- rep(1, ncol(fit$posteriors) - 1)
# names(weights) <-  "Aggregate"
# post <- fit$posteriors[, -1]
# for(i in 1:length(groups)) {
#   group <- groups[[i]]
#   w <- weightList[[1]]
#   sub <- colnames(post) %in% group
#   score <- 1 - apply(post[, sub], 1, function(x) weighted.mean(x, w[sub]))
#   glmfit <- glm(hiv ~ score, family = "binomial")
#   pvals[i] <- summary(glmfit)$coefficients[2, 4] / 2
#   print(roc(hiv ~ score)$auc)
# }
# names(groups) <- paste(names(groups), "pvalue:", round(pvals, 4))

clusterbox <- plot(fit, type = "boxplot", target = infection, groups = groups, ncol = 3,
                   weights = weightList, test = "t-test",
                   one_sided = TRUE)
# save_plot("figures/HVTNisingBoxplotBio.pdf", clusterbox,
#           base_width = 8, base_height = 3)

# Posterior probabilities for nresponses ----------------
infect$hiv <- "non-infected"
infect$hiv[infect$status == 1] <- "infected"
infect$hiv[outcome[, 2] == 1] <- "placebo"

assignments <- fit$assignmentList
names(assignments) <- substr(names(assignments), 1, 9)
assignments <- lapply(unique(names(assignments)), function(x) {
  do.call("rbind", assignments[names(assignments) == x])
})
post <- data.frame(t(sapply(assignments, colMeans)))
post <- cbind(1:nrow(post), post)
subsets <- names(fit$coefficients)

nfunctions <- sapply(subsets, function(x) length(strsplit(x, "+", fixed = TRUE)[[1]]) - 1)
weightList <- list()
weightList$poly <- weights <- nfunctions / choose(5, nfunctions)
weightList$func <- rep(1, length(nfunctions))

resultList <- list()
type <- "poly"
for(w in 1:length(weightList)) {
  weights <- weightList[[w]]
  weightname <- names(weightList)[w]
  resultList[[w]] <- list()
  print(weightname)
  for(i in 1:(length(assignments) + 3)) {
    cat(i, " ")
    if(i <= length(assignments)) {
      samp <- assignments[[i]]
      ptid <- fit$posteriors$ptid[i]
    } else if(i == (length(assignments) + 1)) {
      ptid <- "non-infected"
      samp <- do.call("rbind", assignments[infect[, 3] == "non-infected"])
    } else if(i == (length(assignments) + 2)){
      ptid <- "placebo"
      samp <- do.call("rbind", assignments[infect[, 3] == "placebo"])
    } else {
      ptid <- "infected"
      samp <- do.call("rbind", assignments[infect[, 3] == "infected"])
    }
    #colnames(samp) <- substring(subsets, 1, 3) # for stim groups
    #colnames(samp) <- substring(subsets, 1, 6) # for stim + 8 or 4 groups
    #colnames(samp) <- substring(subsets, 5, 6) # + 8 or 4 groups
    colnames(samp) <- substring(subsets, 5) # for cell-type groups
    #colnames(samp) <- rep("all", ncol(samp)) # for just one plot
    groups <- unique(colnames(samp))
    subjDatList <- list()
    if(type != "breadth") {
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
    } else {
      templist <- list()
      for(j in 1:length(groups)) {
        group <- groups[j]
        subsamp <- samp[, colnames(samp) == group]
        groupSize <- ncol(subsamp)
        weights <- seq(from = 5, to = 1,length.out = groupSize)#(groupSize:1)^0.75
        weights <- weights / sum(weights)
        weights <- cumsum(weights)
        map <- cbind(0:groupSize, c(0, weights))
        values <- map[apply(subsamp, 1, sum), 2]
        templist[[j]] <- values
      }

      templist <- do.call("cbind", templist)
      parents <- substr(groups, 5, 6)
      parentGroups <- c("4+", "8+")
      for(j in 1:length(parentGroups)) {
        parent <- parentGroups[j]
        values <- rowMeans(templist[, parents == parent])
        uniqueVals <- sort(unique(c(values, 0, 1)), decreasing = TRUE)
        props <- sapply(uniqueVals, function(x) mean(values >= x))
        subjDatList[[j]] <- data.frame(ptid = ptid, group = parent,
                                       index = "breadth",
                                       presponses = uniqueVals,
                                       postProb = props)
      }
    }
    resultList[[w]][[i]] <- do.call("rbind", subjDatList)
  }
  cat("\n")
}
responsedat <- do.call("rbind", do.call("rbind", resultList))

forplot <- merge(responsedat, infect, by.x = "ptid", by.y = "ptid",
                 all.x = TRUE)
summarized <- subset(forplot, ptid %in% c("infected", "non-infected", "placebo"))
forplot <- subset(forplot, !(ptid %in% c("infected", "non-infected", "placebo")))
forplot <- forplot[order(forplot$ptid, forplot$group, forplot$presponses), ]

forplot$parent <- substr(forplot$group, 5, 6)
forplot$protein <- substr(forplot$group, 1, 3)
summarized$parent <- substr(summarized$group, 5, 6)
summarized$protein <- substr(summarized$group, 1, 3)
temp <- subset(forplot, index == "poly")
sumtemp <- subset(summarized, index == "poly")
ggplot(temp) + geom_step(aes(x = presponses, y = 1 - postProb, group = ptid,
                                col = factor(hiv)),
                            alpha = 0.42) +
  geom_step(data = sumtemp, aes(x = presponses, y = 1 - postProb,
                                   linetype = ptid), size = 0.7) +
  facet_wrap( ~ group) + theme_bw() +
  xlab("At Least % Responsive Subsets") +
  ylab("Posterior Probabilities") #+

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

intCDF <- summarize(group_by(forplot, ptid, group, index, hiv),
                    measure = integrateCDF(cbind(ptid, group, presponses, postProb)))
intCDF <- do.call("rbind", by(intCDF, list(intCDF$group, intCDF$index),
                              function(x) {
                                measure <- x$measure
                                x$measure <- (measure - mean(measure)) / sd(measure)
                                return(x)
                              }))

par(mfrow = c(1, 2), mar = rep(3, 3))
temp <- subset(intCDF, hiv != "placebo")
aucs <- by(temp, list(temp$group, temp$index), function(x) {
  print(unique(x$group))
  rocfit <- roc(x$hiv == "infected" ~ x$measure)
  auc <- round(rocfit$auc, 3)
  main <- paste(unique(x$group), unique(x$index), "AUC:", auc)
  plot(rocfit, main = main)
  return(auc)
})
n0 <- sum(infect$hiv == "infected")
n1 <- sum(infect$hiv == "non-infected")
pvals <- pwilcox(aucs * n0 * n1, n0, n1, lower.tail = FALSE)

ggplot(subset(intCDF, index != "!!"),
       aes(x = hiv, y = measure, col = hiv)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.1) +
  facet_grid(index ~ group, scales= "free_y") +
  theme_bw() + ylim(-2.5, 2.5)


# Breadth X Functionality Tensor -----------------
library(viridis)
library(ggthemes)
library(dplyr)
library(reshape2)
library(ggplot2)
infect$hiv <- "non-infected"
infect$hiv[infect$status == 1] <- "infected"
infect$hiv[outcome[, 2] == 1] <- "placebo"

post <- fit$posteriors
post <- melt(post, id = "ptid")
names(post)[2:3] <- c("combination", "posterior")
post$subset <- substring(post$combination, 5)
post$stim <- substring(post$combination, 1, 3)
post <- merge(post, infect, all.x = TRUE)
post$hiv <- factor(post$hiv, levels = c("non-infected", "infected", "placebo"))
summ <- summarize(group_by(post, combination, subset, stim, hiv),
                  mposterior = 1 - mean(posterior))
sig <- result$qvals < 0.05
sig <- data.frame(combination = result$subsets, sig = sig)
summ <- merge(summ, sig, all.x = TRUE)
summ$sig[!summ$sig] <- NA
summ$vaccine <- "vaccine"
summ$vaccine[summ$hiv == "placebo"] <- "placebo"
summ$vaccine <- factor(summ$vaccine, level = c("vaccine", "placebo"))

onlysig <- subset(summ, sig)
ggplot(summ, aes(x = stim, y = subset, fill = mposterior)) +
  geom_tile(color = "white", size = 0.1) +
  geom_tile(data = onlysig, color = "red", size = 0.5) +
  scale_fill_viridis(name="posterior") +
  facet_wrap(~ vaccine) +
  theme_tufte(base_family="Helvetica") +
  theme(axis.text.y = element_text(size = 5.5))

sig <- infectresult$infectQvals < 0.1
sig <- data.frame(combination = infectresult$subsets, sig = sig)
summ <- summ[, -which(colnames(summ) == "sig")]
summ <- merge(summ, sig, all.x = TRUE)
summ$sig[!summ$sig] <- NA
temp <- subset(summ, hiv != "placebo")
temp$hiv <- factor(temp$hiv, levels = c("non-infected", "infected"))
onlysig <- subset(temp, sig)
ggplot(temp, aes(x = stim, y = subset, fill = mposterior)) +
  geom_tile(color = "white", size = .1) +
  geom_tile(data = onlysig, aes(x = stim, y = subset),
            col = "red", size = 0.5) +
  scale_fill_viridis(name="posterior") +
  facet_wrap(~ hiv) +
  theme_tufte(base_family="Helvetica") +
  theme(axis.text.y = element_text(size = 5.5))

# Scatter plots for problematic subsets ------------
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

pop <- "8+/107a-154+IFNg+IL2-TNFa+"
pop <- "8+/107a-154+IFNg+IL2+TNFa+"
pop <- "4+/107a+154+IFNg-IL2+TNFa-"
pop <- "4+/107a+154-IFNg-IL2-TNFa+"
pop <- "8+/107a-154+IFNg+IL2-TNFa-"
pop <- "8+/107a-154+IFNg+IL2+TNFa-"
pop <- "8+/107a-154+IFNg+IL2-TNFa-"
pop <- "8+/107a-154+IFNg-IL2-TNFa-"
pop <- "8+/107a+154-IFNg-IL2+TNFa+"
pop <- "8+/107a+154-IFNg-IL2+TNFa+"
pop <- "8+/107a+154-IFNg-IL2+TNFa-"
pop <- "8+/107a+154-IFNg+IL2+TNFa-"
temp <- subset(forplot, population == pop)
ggplot(temp) +
  geom_point(aes(x = log(negprop + 1 / parentcount), y = log(prop + 1 / parentcount),
                 col = 1 - posterior), size = 0.25) +
  facet_wrap( ~ stimGroup, scales = "free", nrow = 2) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_gradientn(colours = rainbow(4)) +
  ggtitle(pop)
sort(temp$prop, decreasing = TRUE)



# Weighting posterior by variances -----
post <- fit$posteriors
post[, 2:ncol(post)] <- 1 - post[2:ncol(post)]
post <- melt(post, id = "ptid")
names(post) <- c("ptid", "subset", "posterior")
post$parent <- substr(post$subset, 5, 6)
post$stim <- substr(post$subset, 1, 3)
post$group <- substr(post$subset, 1, 6)

post$group <- 1
weights <- apply(fit$posteriors[, -1], 2, var)
nfunctions <- sapply(subsets, function(x) length(strsplit(x, "+", fixed = TRUE)[[1]]) - 1)
weights <- weights * nfunctions / choose(5, nfunctions)
weights <- weights / sum(weights)
weights <- data.frame(subset = colnames(fit$posteriors)[-1],
                      weights = weights)
post <- merge(post, weights, all.x = TRUE)
score <- summarize(group_by(post, ptid, group),
                   score = weighted.mean(posterior, weights))
score <- merge(score, infect, all.x = TRUE)

ggplot(score) + geom_boxplot(aes(x = hiv, y = score, col = hiv)) +
  facet_wrap(~ group) + theme_bw()
par(mfrow = c(1, 1))
by(score, score$group, function(x) {
  rocfit <- roc(x$hiv != "placebo" ~ x$score)
  plot(rocfit, main = paste(unique(x$group), round(rocfit$auc, 3)))
})

par(mfrow = c(1, 1))
temp <- subset(score, hiv != "placebo")
by(temp, temp$group, function(x) {
  rocfit <- roc(x$hiv ~ x$score)
  plot(rocfit, main = paste(unique(x$group), round(rocfit$auc, 3)))
})


# Grouping by stim ---------------------------
assignments <- fit$assignmentList
names(assignments) <- substr(names(assignments), 1, 9)
assignments <- lapply(unique(names(assignments)), function(x) {
  do.call("rbind", assignments[names(assignments) == x])
})
sapply(assignments, function(samp) mean(apply(samp, 1, function(x) sum(x == 0) < 40)))
subsets <- names(fit$coefficients)
summarizeBySubset <- function(samp, subsets, threshold = 0) {
  colnames(samp) <- substring(subsets, 5) # for cell-type groups
  #colnames(samp) <- substring(subsets, 1, 6) # for stim groups
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

