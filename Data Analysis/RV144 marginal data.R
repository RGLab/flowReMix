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
data$treatment2 <- data$treatment

preAssignment <- by(data, INDICES = data$ptid, preAssign)
preAssignment <- do.call("rbind", preAssignment)

vaccine <- as.numeric(by(data, data$ptid, function(x) x$vaccine[1] == "VACCINE"))
system.time(fit <- subsetResponseMixtureRcpp(count ~  treatment,
                                         sub.population = factor(data$population),
                                         N = parentcount, id =  ptid,
                                         treatment = treatment2,
                                         data = data,
                                         preAssignment = NULL,
                                         randomAssignProb = 0.0,
                                         weights = NULL,
                                         updateLag = 8, nsamp = 30, maxIter = 16,
                                         isingMethod = "sparse",
                                         covarianceMethod = "dense",
                                         regressionMethod = "betabinom",
                                         centerCovariance = FALSE,
                                         initMHcoef = 3,
                                         dataReplicates = 5,
                                         maxDispersion = 10^2))
#save(fit, file = "dispersed model 2.Robj")
#save(fit, file = "results/binom model.Robj")
#save(fit, file = "results/dispersed model 2 wAssignment.Robj")
#load("results/dispersed model 2 wAssignment.Robj")
#save(fit, file = "dispersed model 3.Robj")
#load(file = "dispersed model 4.Robj")
#load(file = "data analysis/results/marginal independence model od.Robj")
#save(fit, file = "marginal independence model od w covariance.Robj")
save(fit, file = "data analysis/results/marginal model od w ising4.Robj")

require(pROC)
vaccine <- as.vector(by(data, INDICES = data$ptid, FUN = function(x) x$vaccine[1] == "VACCINE"))
posteriors <- fit$posteriors
posteriors$ptid <- as.numeric(as.character(posteriors$ptid))
posteriors <- posteriors[order(posteriors$ptid), -1]
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
  legend('topright', col = c("red", "blue"), lty = 1:2,
         legend = c("FDR", "power"))
}

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
ggplot(forplot) +
  geom_point(aes(x = negprop, y = envprop, col = posterior, shape = vaccine),
             alpha = 0.75) +
  facet_wrap(~ subset, scales = 'free') +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + scale_colour_gradientn(colours=rainbow(4))



