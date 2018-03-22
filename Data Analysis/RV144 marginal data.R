preAssign <- function(dat) {
  subsets <- unique(dat$population)
  nSubsets <- length(subsets)
  preAssign <- numeric(nSubsets)
  prop <- dat$count / dat$parentcount
  for(j in 1:nSubsets) {
    negctrl <- prop[dat$stim == "negctrl" & dat$population == subsets[j]]
    env <- prop[dat$stim == "env" & dat$population == subsets[j]]
    preAssign[j] <- ifelse(env >= negctrl, -1, 0)
  }
  result <- data.frame(id = dat$ptid[1], subset = subsets, assign = preAssign)
  return(result)
}

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
selected_populations = c(1:7)
data <- subset(data, population %in% leaves[selected_populations])
data$population <- factor(data$population)
data <- subset(data, stim != "sebctrl")
data$treatment <- as.numeric(data$stim == "env")
data$ptid <- as.numeric(data$ptid)
data$ptid[data$vaccine == "VACCINE"] <- data$ptid[data$vaccine == "VACCINE"] * 10^4
data$prop <- data$count / data$parentcount
data$population <- as.factor(data$population)
data <- data[order(data$population, data$ptid, data$stim, decreasing = FALSE), ]
data$treatment2 <- data$treatment

control <- flowReMix_control(updateLag = 3, nsamp = 20, initMHcoef = 1,
                             keepEach = 5, isingWprior = FALSE,
                             nPosteriors = 1, centerCovariance = FALSE,
                             maxDispersion = 10^3, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             initMethod = "robust", ncores = 2,
                             markovChainEM = TRUE,
                             seed = 10,
                             preAssignCoefs = 1, sampleNew = FALSE,
                             learningRate = 0.6, keepWeightPercent = 0.9,
                             isingStabilityReps = 0, randStabilityReps = 0,
                             isingInit = -2)

data$stim <- factor(data$stim, levels = c("negctrl", "env"))
# assignmentMat <- do.call("rbind", by(data, data$ptid, preAssign))
system.time(fit <- flowReMix(cbind(count, parentcount - count) ~ stim,
                 subject_id = ptid,
                 cell_type = population,
                 cluster_variable = stim,
                 data = data.frame(data),
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "firth",
                 iterations = 100, parallel = TRUE,
                 cluster_assignment = TRUE, keepSamples = FALSE,
                 verbose = TRUE, control = control, newSampler = FALSE))
# save(fit, file = "Data Analysis/results/RV144 marginals dispersed w all.Robj")
# save(fit, file = "Data Analysis/results/RV144 marginals dispersed wo ising.Robj")
# save(fit, file = "Data Analysis/results/RV144 marginals dispersed wo random.Robj")
# save(fit, file = "Data Analysis/results/RV144 marginals dispersed indepdent.Robj")

system.time(stab <- stabilityGraph(fit, sampleNew = FALSE, reps = 100, cpus = 2))

# S3 methods ----------
plot(stab, fill = roctab$auc, fillRange = c(0.4, 1), plotAll = TRUE,
     label_size = 5, fillPalette = 2, threshold = 1)
plot(fit, type = "scatter", target = vaccine, ncol = 3)
plot(fit, type = "scatter", target = vaccine, ncol = 2)
plot(fit, type = "scatter", target = age, ncol = 2)
plot(fit, type = "scatter", target = vaccine, ncol = 2,
     palette = 2)
plot(fit, type = "scatter", target = vaccine, ncol = 2,
     paletteRange = c(0.5, 1))
plot(fit, type = "scatter", target = vaccine, ncol = 2,
     paletteRange = c(0.5, 1), varname = "flalop")
plot(fit, type = "boxplot", target = vaccine, ncol = 3,
     test = "wilcoxon")
plot(fit, type = "graph", target = vaccine, ncol = 3,
     test = "wilcoxon", fill = roctab$auc)
summary(fit, type = "ROC", target = vaccine)


# Scatter plots -----------------
vaccine <- as.vector(by(data, INDICES = data$ptid, FUN = function(x) x$vaccine[1] == "VACCINE"))
plot(fit, type = "scatter", target = vaccine, ncol = 3)

# ROC table -----------------
roctab <- summary(fit, type = "ROC", target = vaccine)
roctab[order(roctab$auc, decreasing = TRUE), ]

# Plotting Ising --------
plot(fit, type = "graph", threshold = 0, graph = "ising", fill = roctab$auc)
plot(fit, type = "graph", threshold = 0, graph = "randomEffects" ,fill = roctab$auc)

# Plotting Ising --------
plot(fit, type = "graph", threshold = 0, graph = "ising", fill = roctab$auc)
plot(fit, type = "graph", threshold = 0, graph = "randomEffects" ,fill = roctab$auc)

# ROC, FDR and boxplot figures ---------------
plot(fit, type = "ROC", target = vaccine, ncol = 4)
plot(fit, type = "FDR", target = vaccine)
plot(fit, type = "boxplot", target = vaccine,
     test = "wilcoxon", ncol = 4)

# Graphical Models ----------------------
stability <- stabilityGraph(fit, type = "ising", cv = FALSE,
                            reps = 100, cpus = 2)
plot(stability, fill = roctab$auc)

random <- stabilityGraph(fit, type = "randomEffects", cv = FALSE,
                            reps = 50, cpus = 2)
plot(random, fill = roctab$auc, threshold = 0.01)

# ROC, FDR and boxplot figures ---------------
plot(fit, type = "ROC", target = vaccine, ncol = 4)
plot(fit, type = "FDR", target = vaccine)
plot(fit, type = "boxplot", target = vaccine,
     test = "wilcoxon", ncol = 4)

# Graphical Models ----------------------
stability <- stabilityGraph(fit, type = "ising", cv = FALSE,
                            reps = 100, cpus = 2)
plot(stability, fill = roctab$auc)

random <- stabilityGraph(fit, type = "randomEffects", cv = FALSE,
                            reps = 50, cpus = 2)
plot(random, fill = roctab$auc, threshold = 0.01)

# Generarting simulation object ---------
sim <- flowFitToSim(fit)
sim <- generateFlowDataset(sim, 100)
simfit <- fitFlowSim(sim)





