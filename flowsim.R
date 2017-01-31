require(IsingSampler)
require(flowReMix)
load("results/binom model.Robj")
#load("results/dispersed model 2.Robj")
isingmat <- fit$isingCov
randomcov <- fit$covariance
overdispersion <- fit$dispersion

n <- 262
graph <- isingmat
diag(graph) <- 0
thresholds <- diag(isingmat)
assignment <- IsingSampler(n, graph, thresholds)
IsingFit::IsingFit(assignment, AND = FALSE)
rintercept <- mvtnorm::rmvnorm(n, sigma = randomcov)
coefs <- do.call("rbind", fit$coefficients)
batchEffect <- c(0, 0)
subjectlist <- list()
for(i in 1:n) {
  row <- 0
  subjectData <- data.frame(ptid = rep(i, nrow(coefs) * 2))
  controlN <- sample(rv144$parentcount, 1)
  treatmentN <- sample(rv144$parentcount, 1)
  batch <- sample.int(length(batchEffect), 1)
  subjectData$batch <- batch
  for(j in 1:nrow(coefs)) {
    row <- row + 1
    subjectData$treatment[row] <- 0
    subjectData$N[row] <- controlN
    subjectData$eta[row] <- coefs[j, 1] + rintercept[i, j] + batchEffect[batch]
    subjectData$subset[row] <- j
    subjectData$M[row] <- overdispersion[j]

    row <- row + 1
    subjectData$treatment[row] <- 1
    subjectData$N[row] <- treatmentN
    subjectData$eta[row] <- subjectData$eta[row - 1] + coefs[j, 2] * assignment[i, j]
    subjectData$subset[row] <- j
    subjectData$M[row] <- overdispersion[j]
  }
  subjectData$prob <- expit(subjectData$eta)
  subjectData$od <- with(subjectData, rbeta(nrow(subjectData), M * prob, M * (1 - prob)))
  subjectData$count <- rbinom(nrow(coefs) * 2, subjectData$N, subjectData$od)
  subjectlist[[i]] <- subjectData
}
simdata <- do.call("rbind", subjectlist)
simdata$batch <- factor(batch)

system.time(simfit <- subsetResponseMixtureRcpp(count ~  treatment,
                                             sub.population = factor(simdata$subset),
                                             N = N, id =  ptid, treatment = treatment,
                                             data = simdata,
                                             randomAssignProb = 0.0,
                                             rate = 1, updateLag = 3, nsamp = 50, maxIter = 10,
                                             sparseGraph = TRUE,
                                             betaDiserpsion = FALSE,
                                             initMHcoef = 3,
                                             covarianceMethod = c("dense"),
                                             centerCovariance = FALSE))
#save(simfit, file = "results/simfit dispersed")
#save(simfit, file = "results/simfit binomial.Robj")

require(pROC)
posteriors <- 1 - simfit$posteriors[, 2:ncol(fit$posteriors), drop = FALSE]
par(mfrow = c(2, 3), mar = rep(3, 4))
for(i in 1:length(selected_populations)) {
  rocfit <- roc(assignment[, i] ~ posteriors[, i])
  print(plot(rocfit, main = paste(i, "- AUC", round(rocfit$auc, 3))))
}

par(mfrow = c(3, 3), mar = rep(3, 4))
for(i in 1:length(selected_populations)) {
  post <- 1 - posteriors[, i]
  vaccine <- assignment[, i]
  treatment <- assignment[order(post), i]
  uniquePost <- sort(unique(post))
  nominalFDR <- sapply(uniquePost, function(x) mean(post[post <= x]))
  empFDR <- sapply(uniquePost, function(x) 1 - mean(vaccine[post <= x]))
  power <- sapply(uniquePost, function(x) sum(vaccine[post <= x]) / sum(vaccine))
  print(plot(nominalFDR, empFDR, type = "l", xlim = c(0, 1), ylim = c(0, 1), col = "red", main = leaves[selected_populations[i]]))
  lines(nominalFDR, power, col = "blue", lty = 2)
  abline(a = 0, b = 1)
  abline(v = c(0.05, 0.1), h = c(0.8, 0.9), col = "grey")
}

posteriors <- 1 - simfit$posteriors[, 2:ncol(fit$posteriors), drop = FALSE]
forplot <- list()
for(i in 1:length(selected_populations)) {
  post <- 1 - posteriors[, i]
  negprop <- log(simdata$count / simdata$N)[simdata$subset == i & simdata$treatment == 0]
  envprop <- log(simdata$count / simdata$N)[simdata$subset == i & simdata$treatment == 1]
  forplot[[i]] <- data.frame(subset = i,
                             negprop = negprop, envprop = envprop,
                             posterior = 1 - post, vaccine = assignment[, i])
}

forplot <- do.call("rbind", forplot)
require(ggplot2)
ggplot(forplot) +
  geom_point(aes(x = negprop, y = envprop, col = posterior,
                       shape = factor(vaccine))) +
        facet_wrap(~ subset, scales = 'free') +
        geom_abline(slope = 1, intercept = 0) +
        theme_bw() + scale_colour_gradientn(colours=rainbow(4))


# Comparing real data with model data -----------------
selected_populations = c(c(1, 2, 5, 3, 6, 7))
par(mfrow = c(2, 3), mar = rep(3, 4))
for(i in 1:length(selected_populations)) {
  negprop <- log(data$count / data$parentcount)[data$population == leaves[selected_populations[i]] & data$stim == "negctrl"]
  envprop <- log(data$count / data$parentcount)[data$population == leaves[selected_populations[i]] & data$stim == "env"]
  plot(negprop, envprop, pch = ".", col = "red")

  negprop <- log(simdata$count / simdata$N)[simdata$subset == i & simdata$treatment == 0]
  envprop <- log(simdata$count / simdata$N)[simdata$subset == i & simdata$treatment == 1]
  points(negprop, envprop, pch = ".", col = "blue")
  abline(a = 0, b = 1)
}

forplot <- list()
length <- 1
selected_populations <- c(1, 2, 5, 3, 6, 7)
for(i in 1:length(selected_populations)) {
  negprop <- log(data$count / data$parentcount)[data$population == leaves[selected_populations[i]] & data$stim == "negctrl"]
  envprop <- log(data$count / data$parentcount)[data$population == leaves[selected_populations[i]] & data$stim == "env"]
  forplot[[length]] <- data.frame(negprop = negprop, envprop = envprop,
                                  subpop = leaves[selected_populations[i]],
                                  data = "rv144")
  length <- length + 1

  negprop <- log(simdata$count / simdata$N)[simdata$subset == i & simdata$treatment == 0]
  envprop <- log(simdata$count / simdata$N)[simdata$subset == i & simdata$treatment == 1]
  forplot[[length]] <- data.frame(negprop = negprop, envprop = envprop,
                                  subpop = leaves[selected_populations[i]],
                                  data = "sim")
  length <- length + 1
}

forplot <- do.call("rbind", forplot)
ggplot(forplot, aes(x = negprop, y = envprop, col = data)) +
  geom_density_2d(aes(linetype = data)) +
  theme_bw() + facet_wrap(~ subpop) + geom_abline(intercept = 0, slope = 1)
