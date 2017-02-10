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
selected_populations = c(1:4, 5:7)
data <- subset(data, population %in% leaves[selected_populations])
data$population <- factor(data$population)
data <- subset(data, stim != "sebctrl")
data$treatment <- as.numeric(data$stim == "env")
data$ptid <- as.numeric(data$ptid)
data$ptid[data$vaccine == "VACCINE"] <- data$ptid[data$vaccine == "VACCINE"] * 10^4
data$prop <- data$count / data$parentcount
data$population <- as.factor(data$population)
data <- data[order(data$population, data$ptid, data$stim, decreasing = FALSE), ]

system.time(fit <- flowRegressionMixture(count ~  treatment,
                      sub.population = factor(data$population),
                      N = parentcount, id =  ptid,
                      data = data,
                      treatment = treatment,
                      weights = NULL,
                      rate = 1, updateLag = 5,
                      nsamp = 200,
                      centerCovariance = FALSE,
                      maxIter = 20, tol = 1e-03))

# Facet Wrap Plot! ---------------------------
posteriors <- fit$posteriors[, 3]
populations <- unique(data$population)
plotList <- lapply(1:length(selected_populations), function(x) x)
for(i in 1:length(populations)) {
  tempdat <- subset(data, population == populations[i])
  prop <- log(tempdat$count / tempdat$parentcount)
  stimprop <- prop[tempdat$stim == "env"]
  ctrlprop <- prop[tempdat$stim != "env"]
  vaccine <- (tempdat$vaccine == "VACCINE")[tempdat$stim == "env"]
  plotList[[i]] <- data.frame(vaccine = vaccine, stimprop = stimprop,
                              ctrlprop = ctrlprop, population = populations[i],
                              posterior = posteriors)
}

forPlot <- do.call("rbind", plotList)
require(ggplot2)
ggplot(forPlot) +
  geom_point(aes(x = ctrlprop, y = stimprop, col = posterior, shape = !vaccine), fill = "white") +
  theme_bw() + geom_abline(intercept = 0, slope = 1) +
  scale_colour_gradientn(colours=rainbow(4)) +
  facet_wrap(~ population, scales = "free", ncol = 4)


par(mfrow = c(1, 2), mar = rep(3, 4))
sprobs <- data.frame(posteriors, vaccine)
sprobs <- sprobs[order(posteriors, decreasing = TRUE), ]
sprobs$nominalFDR <- cummean(1 - sprobs$posteriors)
sprobs$empFDR <- cummean(1 - sprobs$vaccine)
sprobs$power <- cumsum(sprobs$vaccine) / sum(sprobs$vaccine)
uniqueNominal <- unique(sprobs$nominalFDR)
empFDR <- sapply(uniqueNominal, function(x) sprobs$empFDR[max(which(sprobs$nominalFDR == x))])
power <- sapply(uniqueNominal, function(x) sprobs$power[max(which(sprobs$nominalFDR == x))])
lim <- max(c(empFDR, uniqueNominal))
plot(uniqueNominal, empFDR, type = "l", xlim = c(0, lim), ylim = c(0, 1), col = "red",
     xlab = "nominal FDR", ylab = "Empirical FDR / Power")
abline(a = 0, b = 1)
lines(uniqueNominal, power, col = "blue", lty = 2)
legend("topright", col = c("red", "blue"), lty = 1:2, legend = c("FDR", "Power"))
abline(v = c(.01, .05, .1), h = c(.75, .9, .95), col = "grey")
rocfit <- roc(vaccine ~ posteriors)
print(plot(rocfit, main = round(rocfit$auc, 3)))

covariance <- fit$covariance
require(xtable)
covTable <- xtable(covariance, digits = 2)
names(covTable) <- 1:nlevels(data$population)
rownames(covTable) <- leaves[c(1:3,5:7)]
corMat <- cov2cor(covariance)
corTable <- xtable(corMat, digits = 2)
names(corTable) <- 1:ncol(corTable)
rownames(corTable) <- leaves[c(1:3,5:7)]

# Malaria dataset
data <- malaria
leaves <- unique(data$population)
leaves <- leaves[!(leaves %in% unique(data$parent))]
leaf <- leaves[2:4]
data <- subset(data, population %in% leaf)
data$visitno <- as.numeric(factor(data$visitno))
data$visitno <- factor(data$visitno - min(data$visitno))
data$treatment <- data$visitno
fit <- flowRegressionMixture(count ~  (treatment * stim + experiment),
                             sub.population = factor(data$population),
                             N = parentcount, id =  ptid,
                             data = data,
                             treatment = treatment,
                             weights = NULL,
                             rate = 1, updateLag = 5,
                             nsamp = 200,
                             centerCovariance = TRUE,
                             maxIter = 30, tol = 1e-03)
sapply(fit$coefficients, function(x) x)






