library(flowReMix)
library(pROC)
library(MIMOSA)
cummean <- function(x) cumsum(x) / 1:length(x)
data(rv144)
omit <- paste("P", c(1001, 1013, 1019, 1023, 1031, 1034, 1039, 1045,
                     1060, 1095, 1099, 1100, 1109, 1177, 1180, 1187,
                     1201, 1215, 1216, 1224, 1227, 1232, 1242, 1284),
              sep = "")
par(mfrow = c(1, 1), mar = rep(4, 4))
data <- rv144
data <- subset(data, !(ptid %in% omit))
leaves <- unique(data$population)
selected_populations = c(c(1, 2, 5, 3, 6, 7))
data <- subset(data, population %in% leaves[selected_populations])
data$population <- factor(data$population)
data <- subset(data, stim != "sebctrl")
data$treatment <- as.numeric(data$stim == "env")
data$ptid <- as.numeric(data$ptid)
data$ptid[data$vaccine == "VACCINE"] <- data$ptid[data$vaccine == "VACCINE"] * 10^4
data$prop <- data$count / data$parentcount
data <- data[order(data$population, data$ptid, data$stim, decreasing = FALSE), ]

# Analysis with MIMOSA
library(MIMOSA)
count <- data$count
negcount <- data$parentcount - data$count
antigen <- data$stim
cytokine <- data$population
ptid <- data$ptid
mimosaData <- data.frame(count, negcount, antigen, cytokine, ptid)
mimosaObject <- ConstructMIMOSAExpressionSet(mimosaData,
                                             reference = antigen %in% "negctrl",
                                             measure.columns = c("count", "negcount"),
                                             other.annotations = c("cytokine", "antigen", "ptid"),
                                             .variables = .(cytokine, ptid),
                                             default.cast.formula = component ~ ptid + antigen + cytokine)
mimosafit <- MIMOSA(negcount + count ~ ptid + antigen | cytokine,
                    data = mimosaObject, method = "EM",subset = RefTreat %in% "Treatment",
                    ref = RefTreat %in% "Reference")

load(file = "data analysis/results/marginal independence model od.Robj")

# AUC comparison
posteriors <- fit$posteriors[order(as.numeric(as.character(fit$posteriors$ptid))), -1]
vaccine <- as.numeric(by(data, data$ptid, function(x) x$vaccine[1] == "VACCINE"))
par(mfrow = c(2, 3), mar = rep(1, 4))
for(j in 1:length(mimosafit)) {
  subset <- names(mimosafit)[j]
  mixIndex <- which(names(posteriors) == subset)

  mimosaPvals <- getZ(mimosafit[[j]])[, 2]
  mixposteriors <- posteriors[, mixIndex]

  mimosaroc <- roc(vaccine ~ mimosaPvals)
  mixroc <- roc(vaccine ~ mixposteriors)
  plot(mimosaroc, col = "red", lty = 2,
       main = paste(subset, "mimosa:", round(mimosaroc$auc, 3), " - reg:", round(mixroc$auc, 3)))
  lines(mixroc, col = "black", lty = 4)
  legend("bottomright", col = c("red", "black"), lty = c(2, 4),
         legend = c("MIMOSA", "Regression"))
}

# FDR comparison
par(mfrow = c(2, 3), rep(4, 4))
for(j in 1:length(mimosafit)) {
  subset <- names(mimosafit)[j]
  mixIndex <- which(names(posteriors) == subset)

  mimosaQvals <-  fdr(mimosafit[[j]])
  mixposteriors <- posteriors[, mixIndex]

  post <-  posteriors[, j]
  treatment <- vaccine[order(post)]
  uniquePost <- sort(unique(post))
  nominalFDR <- sapply(uniquePost, function(x) mean(post[post <= x]))
  empFDR <- sapply(uniquePost, function(x) 1 - mean(vaccine[post <= x]))
  plot(nominalFDR, empFDR, type = "l", xlim = c(0, 1), ylim = c(0, 1), col = "black",
       main = subset)
  abline(a = 0, b = 1, col = "grey")

  mimosaFDR <- sapply(uniquePost, function(x) 1 - sum(vaccine[mimosaQvals <= x]) / sum(mimosaQvals <= x))
  lines(nominalFDR, mimosaFDR, col = "red")
  legend("topright", col = c("red", "black"), lty = 1,
         legend = c("MIMOSA", "Regression"))
}


