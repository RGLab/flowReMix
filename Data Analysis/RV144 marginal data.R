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

control <- flowReMix_control(updateLag = 10, nsamp = 50, initMHcoef = 1,
                             nPosteriors = 1, centerCovariance = TRUE,
                             maxDispersion = 10^3 / 2, minDispersion = 10^6,
                             randomAssignProb = 0.2, intSampSize = 50,
                             initMethod = "binom", ncores = NULL)

system.time(fit <- flowReMix(cbind(count, parentcount - count) ~ treatment + age + gender,
                 subject_id = ptid,
                 cell_type = population,
                 cluster_variable = treatment,
                 data = data,
                 covariance = "sparse",
                 ising_model = "raIsing",
                 regression_method = "betabinom",
                 iterations = 5, parallel = TRUE,
                 verbose = TRUE, control = control))
#save(fit, file = "Data Analysis/results/RV144 marginals dispersed model new 2.Robj")

## ROC ------------------
require(pROC)
vaccine <- as.vector(by(data, INDICES = data$ptid, FUN = function(x) x$vaccine[1] == "VACCINE"))
ids <- factor(unique(data$ptid), levels = fit$posteriors$ptid)
posteriors <- fit$posteriors
vaccine <- vaccine[order(ids)]
par(mfrow = c(3, 3), mar = rep(3, 4))
aucs <- numeric(ncol(posteriors) - 1)
for(i in 2:ncol(posteriors)) {
  rocfit <- roc(vaccine ~ posteriors[, i])
  print(plot(rocfit, main = paste(colnames(posteriors)[i], "- AUC", round(rocfit$auc, 3))))
  aucs[i - 1] <- rocfit$auc
}
n0 <- sum(vaccine == 0)
n1 <- sum(vaccine == 1)
pvals <- pwilcox(aucs * n1 * n0, n0, n1, lower.tail = FALSE)
qvals <- p.adjust(pvals, method = "BY")

# Scatter plots --------------------------
forplot <- list()
vaccine <- as.vector(by(data, INDICES = data$ptid, FUN = function(x) x$vaccine[1] == "VACCINE"))
ids <- factor(unique(data$ptid), levels = fit$posteriors$ptid)
posteriors <- fit$posteriors
posteriors$ptid <- as.numeric(as.character(fit$posteriors$ptid))
posteriors <- posteriors[order(posteriors$ptid), ]
for(i in 1:length(selected_populations)) {
  post <- posteriors[, i + 1]
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

# Infection status ---------------
data("rv144_correlates_data")
vaccine <- as.vector(by(data, INDICES = data$ptid, FUN = function(x) x$vaccine[1] == "VACCINE"))
correlates <- rv144_correlates_data
correlates <- correlates[order(correlates$PTID), ]
infection <- correlates$infect.y
posteriors <- fit$posteriors
posteriors$ptid <- as.numeric(as.character(fit$posteriors$ptid))
posteriors <- posteriors[order(posteriors$ptid), ]
par(mfrow = c(4, 4), mar = rep(3, 4))
aucs <- numeric(ncol(posteriors) - 1)
par(mfrow = c(3, 3), mar = rep(3, 4))
for(i in 2:ncol(posteriors)) {
  rocfit <- roc(infection[vaccine] ~ posteriors[vaccine, i])
  print(plot(rocfit, main = paste(colnames(posteriors)[i], "- AUC", round(rocfit$auc, 3))))
  aucs[i - 1] <- rocfit$auc
}

n0 <- sum(infect == "infected")
n1 <- sum(infect == "non-infected")
pvals <- pwilcox(aucs * n1 * n0, n0, n1, lower.tail = FALSE)


# Testing isotonic regression -----------
assignments <- fit$assignmentList
names(assignments) <- sapply(names(assignments), function(x) strsplit(x, "%%%")[[1]][[1]])
assignments <- lapply(unique(names(assignments)), function(x) {
  do.call("rbind", assignments[names(assignments) == x])
})

samp <- do.call("rbind", assignments)
samp <- rowSums(samp)
mle <- table(samp)
mle <- mle / sum(mle)
pone <- sort(mle, decreasing = TRUE)
pone <- pone / sum(pone)
plot(diff(log(pone)))




