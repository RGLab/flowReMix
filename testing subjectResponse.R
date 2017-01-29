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
data$population <- as.factor(data$population)
data <- data[order(data$population, data$ptid, data$stim, decreasing = FALSE), ]

vaccine <- as.numeric(by(data, data$ptid, function(x) x$vaccine[1] == "VACCINE"))
system.time(fit <- subsetResponseMixtureRcpp(count ~  treatment,
                                         sub.population = factor(data$population),
                                         N = parentcount, id =  ptid, treatment = treatment,
                                         data = data,
                                         randomAssignProb = 0.0,
                                         weights = NULL,
                                         rate = 1, updateLag = 3, nsamp = 50, maxIter = 12,
                                         sparseGraph = TRUE,
                                         covarianceMethod = c("dense"),
                                         centerCovariance = FALSE))

# system.time(fit <- subsetResponseMixtureNested(count ~  treatment,
#                                              sub.population = factor(data$population),
#                                              N = parentcount, id =  ptid,
#                                              data = data,
#                                              treatment = treatment,
#                                              weights = NULL,
#                                              rate = 1, updateLag = 3,
#                                              nsamp = 40,
#                                              centerCovariance = FALSE,
#                                              maxIter = 8, tol = 1e-03))


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
  abline(v = c(0.05, 0.1), h = c(0.8, 0.9), col = "grey")
}

posteriors <- fit$posteriors[, 2:ncol(fit$posteriors), drop = FALSE]
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
print(ggplot(forplot) +
  geom_point(aes(x = negprop, y = envprop, col = posterior, shape = vaccine)) +
  facet_wrap(~ subset, scales = 'free') +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + scale_colour_gradientn(colours=rainbow(4)))



# Booleans dataset --------------------------------------------
data("rv144_booleans")
bySubset <- by(data.frame(booleans$stim, booleans$nonstim), booleans$Subset, function(x) x)
largerThanThershold <- sapply(bySubset, function(x) colSums(x >5))

require(reshape2)
booldata <- melt(booleans, c("PTID", "Subset"))
names(booldata)[3:4] <- c("stim", "count")

forParentcount <- rv144
forParentcount <- as.data.frame(forParentcount)
forParentcount <- subset(forParentcount,
                         parent == "4+" & stim == "env")
forParentcount <- forParentcount[, c(2, 5, 11)]
forParentcount <- unique(forParentcount)

booldata <- merge(booldata, forParentcount, by.x = "PTID", by.y = "ptid",
                  all.x = TRUE, all.y = FALSE)
booldata <- subset(booldata, Subset != "!TNFa&!IFNg&!IL4&!IL2&!CD154&!IL17a")
booldata$treatment <- as.numeric(booldata$stim == "stim")
uniquepop <- unique(booldata$Subset)
booldata <- with(booldata, booldata[order(Subset, PTID, stim, decreasing = FALSE), ])
booldata <- subset(booldata, !is.na(Subset))
allsubset <- booldata
booldata <- with(booldata, booldata[order(Subset, PTID, stim, decreasing = FALSE), ])

vaccine <- as.numeric(by(booldata, booldata$PTID, function(x) x$vaccine[1] == "VACCINE"))
require(pROC)
system.time(fit <- subsetResponseMixtureRcpp(count ~  treatment,
                                             sub.population = factor(booldata$Subset),
                                             N = parentcount, id =  PTID,
                                             data = booldata,
                                             treatment = treatment,
                                             randomAssignProb = 0,
                                             weights = NULL,
                                             rate = 1, updateLag = 3,
                                             nsamp = 50, maxIter = 15,
                                             covarianceMethod = "sparse",
                                             sparseGraph = TRUE,
                                             centerCovariance = FALSE))
subsets <- unique(booldata$Subset)

require(pROC)
posteriors <- fit$posteriors[, 2:ncol(fit$posteriors), drop = FALSE]
par(mfrow = c(4, 6), mar = rep(1, 4))
auc <- numeric(length(subsets))
for(i in 1:length(subsets)) {
  try(rocfit <- roc(!vaccine ~ posteriors[, i]))
  try(rocfit <- roc(vaccine ~ posteriors[, i]))
  auc[i] <- rocfit$auc
  print(plot(rocfit, main = paste(i, "- AUC", round(rocfit$auc, 3))))
}


par(mfrow = c(4, 6), mar = rep(2, 4))
for(i in 1:length(subsets)) {
  post <- posteriors[, i]
  treatment <- vaccine[order(post)]
  uniquePost <- sort(unique(post))
  nominalFDR <- sapply(uniquePost, function(x) mean(post[post <= x]))
  empFDR <- sapply(uniquePost, function(x) 1 - mean(vaccine[post <= x]))
  power <- sapply(uniquePost, function(x) sum(vaccine[post <= x]) / sum(vaccine))
  print(plot(nominalFDR, empFDR, type = "l", xlim = c(0, 1), ylim = c(0, 1), col = "red",
             main = paste(i, "- AUC ", round(auc[i], 2))))
  lines(nominalFDR, power, col = "blue", lty = 2)
  abline(a = 0, b = 1)
  abline(v = c(0.05, 0.1), h = c(0.8, 0.9), col = "grey")
}


forplot <- list()
for(i in 1:length(subsets)) {
  post <- posteriors[, i]
  negprop <- log(booldata$count / booldata$parentcount)[booldata$Subset == subsets[i] & booldata$stim == "nonstim"]
  envprop <- log(booldata$count / booldata$parentcount)[booldata$Subset == subsets[i] & booldata$stim == "stim"]
  forplot[[i]] <- data.frame(subset = subsets[i],
                             negprop = negprop, envprop = envprop,
                             posterior = 1 - post, vaccine = vaccine)
}

forplot <- do.call("rbind", forplot[c(1:21, 23)])
require(ggplot2)
print(ggplot(forplot) +
        geom_point(aes(x = negprop, y = envprop, col = posterior, shape = vaccine == 1)) +
        facet_wrap(~ subset, scales = 'free', ncol = 6) +
        geom_abline(slope = 1, intercept = 0) +
        theme_bw() + scale_colour_gradientn(colours=rainbow(4)))



library(gridExtra)
table <- data.frame(index = 1:length(subsets), subset = subsets,
                    responseProb = round(fit$levelProbs,2),
                    AUC = round(auc, 2))
pdf("figures/cell subset table B.pdf", height=8, width=6)
grid.table(table, rows = NULL)
dev.off()



