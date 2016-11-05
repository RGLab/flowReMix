# Mixed Effects ----------
data <- rv144
leaves <- unique(data$population)
data <- subset(data, stim != "sebctrl")
data$treatment <- as.numeric(data$stim == "env")
leaf <- leaves[c(1, 7)]
data <- subset(data, population %in% leaf)
data$ptid <- as.numeric(data$ptid)
data$ptid[data$vaccine == "VACCINE"] <- data$ptid[data$vaccine == "VACCINE"] * 10^4
data$prop <- data$count / data$parentcount
data$population <- as.factor(data$population)
data$stim[data$stim == "env"] <- "env"
data$stimB <- data$stim == "env"
data <- data[order(data$ptid, decreasing = FALSE), ]

fit <- lme4::glmer(cbind(count, parentcount - count) ~ age + gender + stim*vaccine + (1|ptid),
             family = binomial,
             data = data)
predict1 <- predict(fit)
predict0 <- predict(fit, re.form = NA)


system.time(fit <- lme4::glmer(cbind(count, parentcount - count) ~ population/(age + gender + stim*vaccine)
                         + (0 + population|ptid),
             family = binomial,
             data = data,
             nAGQ = 1))
summary(fit)
library(brms)
system.time(fit <- brm(count ~ population/(age + gender + stim*vaccine)
                         + (1|population:ptid),
                         family = 'binomial',
                         data = data))
library(glmmstan)
system.time(fit <- glmmstan(cbind(count, parentcount - count) ~ population/(age + gender + stim*vaccine)
                       + (0 + population|ptid),
                       family = "beta-binomial",
                       data = data))

library(MASS)
system.time(fit <- glmmPQL(fixed = cbind(count, parentcount - count) ~ population/(age + gender + stim*vaccine),
                            random = ~ population | ptid,
                            family = binomial,
                            data = data,
                           niter = 100))


# Mixture model -----------------------
data <- rv144
leaves <- unique(data$population)
data <- subset(data, stim != "sebctrl")
data$treatment <- as.numeric(data$stim == "env")
leaf <- leaves[c(2)]
data <- subset(data, population %in% leaf)
data$ptid <- as.numeric(data$ptid)
data$ptid[data$vaccine == "VACCINE"] <- data$ptid[data$vaccine == "VACCINE"] * 10^4
data$prop <- data$count / data$parentcount
data$population <- as.factor(data$population)
data <- data[order(data$ptid, decreasing = FALSE), ]
#data$treatment <- as.numeric(data$vaccine == "VACCINE")

mixedfit <- glmmMixture(count ~ (age + gender + treatment),
                        N = parentcount,
                        sub.population = population,
                        id = ptid,
                        treatment = treatment,
                        data = data,
                        tol = 0.002,
                        maxiter = 20,
                        nAGQ = 1)

cummean <- function(x) cumsum(x) / 1:length(x)
sprobs <- data.frame(nullprob = mixedfit$subject.probs[, 1], altprob = mixedfit$subject.probs[, 2])
rowIndex <- seq(from = nrow(data) / 286, to = nrow(data), by = nrow(data) / 286)
sprobs$vaccine <- data$vaccine[rowIndex]
sprobs$id <- data$ptid[rowIndex]
sprobs$nullprop <- log(data$count/data$parentcount + 10^-5)[data$stim == "negctrl" & data$population == leaf[1]]
sprobs$altprop <- log(data$count/data$parentcount + 10^-5)[data$stim == "env" & data$population == leaf[1]]
sprobs$nonrespNullProp <- log(mixedfit$mu$nullMu[mixedfit$mu$waves == 1] + 10^-5)
sprobs$nonrespAltProp <- log(mixedfit$mu$nullMu[mixedfit$mu$waves == 2] + 10^-5)
sprobs$respNullProp <- log(mixedfit$mu$altMu[mixedfit$mu$waves == 1] + 10^-5)
sprobs$respAltProp <- log(mixedfit$mu$altMu[mixedfit$mu$waves == 2] + 10^-5)
sprobs <- sprobs[order(sprobs$nullprob), ]
sprobs$empFDR <- cummean(sprobs$vaccine == "PLACEBO")
sprobs$nominalFDR <- cummean(sprobs$nullprob)
sprobs$power <- cumsum(sprobs$vaccine == "VACCINE") / sum(sprobs$vaccine == "VACCINE")


par(mfrow = c(1 , 2), mar = rep(4, 4))
uniqueNominal <- unique(sprobs$nominalFDR)
empFDR <- sapply(uniqueNominal, function(x) sprobs$empFDR[max(which(sprobs$nominalFDR == x))])
power <- sapply(uniqueNominal, function(x) sprobs$power[max(which(sprobs$nominalFDR == x))])
lim <- max(c(empFDR, uniqueNominal))
plot(uniqueNominal, empFDR, type = "l", xlim = c(0, lim), ylim = c(0, 1), col = "red",
     xlab = "nominal FDR", ylab = "Empirical FDR / Power")
abline(a = 0, b = 1)
lines(uniqueNominal, power, col = "blue", lty = 2)
legend("topright", col = c("red", "blue"), lty = 1:2, legend = c("FDR", "Power"))

require(pROC)
rocfit <- pROC::roc(sprobs$vaccine ~ sprobs$altprob)
plot(rocfit, main = paste(leaf[1], "    AUC:", round(rocfit$auc, 3)))

sprobs$vaccine <- sprobs$vaccine == "VACCINE"
require(ggplot2)
ggplot(sprobs) +
  geom_point(aes(x = nullprop, y = altprop, col = altprob, shape = !vaccine, size = !vaccine), fill = "white") +
  theme_bw() + geom_abline(intercept = 0, slope = 1) +
  scale_colour_gradientn(colours=rainbow(4)) +
  stat_smooth(aes(x = nullprop, y = nonrespAltProp), col = "red") +
  stat_smooth(aes(x = nullprop, y = respAltProp), col = "purple")

# Comparison with MIMOSA ----------------------
mimosaList <- list()
glmList <- list()
par(mfrow = c(2, 4), mar = rep(3, 4))
set.seed(500)
for(i in 1:7) {
  data <- rv144
  leaves <- unique(data$population)
  data <- subset(data, stim != "sebctrl")
  data$treatment <- as.numeric(data$stim == "env")
  leaf <- leaves[c(i)]
  data <- subset(data, population %in% leaf)
  data$ptid <- as.numeric(data$ptid)
  data$ptid[data$vaccine == "VACCINE"] <- data$ptid[data$vaccine == "VACCINE"] * 10^4
  data$prop <- data$count / data$parentcount
  data$population <- as.factor(data$population)
  data$stim[data$stim == "env"] <- "env"
  data$stimB <- data$stim == "env"
  data <- data[order(data$ptid, decreasing = FALSE), ]

  count <- data$count
  negcount <- data$parentcount - count
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
  mimosafit <- MIMOSA(negcount + count ~ ptid + antigen | cytokine, data= mimosaObject, method = "EM",subset=RefTreat%in%"Treatment",ref=RefTreat%in%"Reference")
  mimosaList[[i]] <- mimosafit

  mimosaQvals <- fdr(mimosafit[[1]])

  mixedfit <- glmmMixture(count ~ (treatment),
                          N = parentcount,
                          sub.population = population,
                          id = ptid,
                          treatment = treatment,
                          data = data,
                          tol = 0.002,
                          maxiter = 20,
                          nAGQ = 1)
  glmList[[i]] <- mixedfit

  glmRoc <- roc(data$vaccine[data$stim == "env"] ~ mixedfit$subject.probs[, 2])
  mimosaRoc <- roc(data$vaccine[data$stim == "env"] ~ mimosaQvals)
  plot(glmRoc, col = "red", main = leaf)
  lines(mimosaRoc, col = "blue")
  legend("bottomright", col = c("red", "blue"), lty = 1,
         legend = c(paste("glmer", round(glmRoc$auc, 3)), paste("mimosa", round(mimosaRoc$auc, 3))))
}

# Plotting FDR
par(mfrow = c(2, 4), mar = rep(4, 4))
for(i in 1:7) {
  mimosafit <- mimosaList[[i]]
  mimosaQvals <- fdr(mimosafit[[1]])
  empmimosaFDR <- sapply(mimosaQvals, function(x) mean(!vaccine[mimosaQvals <= x]))
  mimosapower <- sapply(mimosaQvals, function(x) sum(vaccine[mimosaQvals <= x])) / sum(vaccine)

  vaccine <- data$vaccine[data$stim == "env"] == "VACCINE"
  mixedfit <- glmList[[i]]
  mixqvals <- mixedfit$subject.probs[, 1]
  mixqvals <- sapply(mixqvals, function(x) mean(mixqvals[mixqvals <= x]))
  empmixFDR <- sapply(mixqvals, function(x) mean(!vaccine[mixqvals <= x]))
  mixpower <- sapply(mixqvals, function(x) sum(vaccine[mixqvals <= x])) / sum(vaccine)
  order <- order(mixqvals)
  plot(mixqvals[order], empmixFDR[order], type = "l", col = "red",
       xlim = c(0, 1), ylim = c(0, 1), main = leaves[[i]],
       xlab = "nominal FDR", ylab = "empricial FDR")
  #lines(mixqvals[order], mixpower[order], col = "red", lty = 2)
  abline(a = 0, b = 1)
  order <- order(mimosaQvals)
  lines(mimosaQvals[order], empmimosaFDR[order], col = "blue")
  #lines(mimosaQvals[order], mimosapower[order], col = "blue", lty = 2)
  legend("topright", legend = c("glmer", "MIMOSA"), col = c("red", "blue"),
         lty = 1)
}


# Mapping the Correlations ------------------------
allPairs <- t(combn(c(1:7), 2))
data <- rv144
leaves <- unique(data$population)
data <- subset(data, stim != "sebctrl")
data$treatment <- as.numeric(data$stim == "env")
data$ptid <- as.numeric(data$ptid)
data$ptid[data$vaccine == "VACCINE"] <- data$ptid[data$vaccine == "VACCINE"] * 10^4
data$prop <- data$count / data$parentcount
data$population <- as.factor(data$population)
data$stim[data$stim == "env"] <- "env"
data$stimB <- data$stim == "env"
data <- data[order(data$ptid, decreasing = FALSE), ]
covMat <- matrix(ncol = 7, nrow = 7)
diag(covMat) <- 1
for(i in 1:nrow(allPairs)) {
  pair <- allPairs[i, ]
  fit <- lme4::glmer(cbind(count, parentcount - count) ~ population/(age + gender + stim*vaccine)
                                 + (0 + population|ptid),
                                 family = binomial,
                                 data = subset(data, population %in% leaves[pair]))
  summ <- summary(fit)
  pair[pair == 15] <- 8
  covMat[pair, pair] <- summ$varcor$ptid[1:2, 1:2]
  print(i / nrow(allPairs))
  print(round(covMat, 3))
}
#save(covMat, file = "figures/covMat.Robj")
require(xtable)
covTable <- xtable(covMat, digits = 2)
rownames(covTable) <- leaves[1:7]
#names(covTable) <- leaves[1:7]
corMat <- cov2cor(covMat)
corTable <- xtable(corMat, digits = 2)
rownames(corTable) <- leaves[1:7]







