require(IsingSampler)
isingmat <- fit$isingCov
randomcov <- fit$covariance

n <- 262
graph <- isingmat
diag(graph) <- 0
thresholds <- diag(isingmat)
assignment <- IsingSampler(n, graph, thresholds)
rintercept <- mvtnorm::rmvnorm(n, sigma = randomcov)
coefs <- do.call("rbind", fit$coefficients)
batchEffect <- c(0, 0)
subjectlist <- list()
M <- 5000
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

    row <- row + 1
    subjectData$treatment[row] <- 1
    subjectData$N[row] <- treatmentN
    subjectData$eta[row] <- subjectData$eta[row - 1] + coefs[j, 2] * assignment[i, j]
    subjectData$subset[row] <- j
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
                                             rate = 1, updateLag = 3, nsamp = 40, maxIter = 8,
                                             sparseGraph = TRUE,
                                             covarianceMethod = c("dense"),
                                             centerCovariance = FALSE))

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


