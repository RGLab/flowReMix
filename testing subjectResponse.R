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
selected_populations = c(c(1:2))
data <- subset(data, population %in% leaves[selected_populations])
data$population <- factor(data$population)
data <- subset(data, stim != "sebctrl")
data$treatment <- as.numeric(data$stim == "env")
data$ptid <- as.numeric(data$ptid)
data$ptid[data$vaccine == "VACCINE"] <- data$ptid[data$vaccine == "VACCINE"] * 10^4
data$prop <- data$count / data$parentcount
data$population <- as.factor(data$population)
data <- data[order(data$population, data$ptid, data$stim, decreasing = FALSE), ]

system.time(fit <- subsetResponseMixture(count ~  treatment,
                                         sub.population = factor(data$population),
                                         N = parentcount, id =  ptid,
                                         data = data,
                                         treatment = treatment,
                                         weights = NULL,
                                         rate = 1, updateLag = 5,
                                         nsamp = 30,
                                         centerCovariance = FALSE,
                                         maxIter = 20, tol = 1e-03))

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
}

forplot <- list()
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
  theme_bw())



