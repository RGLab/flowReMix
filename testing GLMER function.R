library(flowReMix)
cummean <- function(x) cumsum(x) / 1:length(x)
data(rv144)
#set.seed(502)
par(mfrow = c(1, 1), mar = rep(4, 4))
data <- rv144
leaves <- unique(data$population)
selected_populations = c(1:3, 5:7)
data <- subset(data, population %in% leaves[selected_populations])
data$population=factor(data$population)
data <- subset(data, stim != "sebctrl")
data$treatment <- as.numeric(data$stim == "env")
data$ptid <- as.numeric(data$ptid)
data$ptid[data$vaccine == "VACCINE"] <- data$ptid[data$vaccine == "VACCINE"] * 10^4
data$prop <- data$count / data$parentcount
data$population <- as.factor(data$population)
data <- data[order(data$population, data$ptid, data$stim, decreasing = FALSE), ]

fit <- flowRegressionMixture(count ~  treatment,
                      sub.population = factor(data$population),
                      N = parentcount, id =  ptid,
                      data = data,
                      treatment = treatment,
                      weights = NULL,
                      rate = 1,
                      nsamp = 100,
                      maxIter = 50, tol = 1e-03)

# Facet Wrap Plot! ---------------------------
posteriors <- fit$posteriors[, 3]
populations <- unique(data$population)
plotList <- lapply(1:length(mixtureFitList), function(x) x)
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
  facet_wrap(~ population, scales = "free", ncol = 2)






