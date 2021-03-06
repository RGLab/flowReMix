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

require(pROC)
pdf("booleanOutputPDF.pdf", paper = "letter")
for(i in 1:length(uniquepop)) {
  booldata <- allsubset# subset(allsubset, Subset %in% uniquepop[c(i)])
  require(flowReMix)
  fit <- flowRegressionMixture(count ~  treatment,
                                 sub.population = factor(booldata$Subset),
                                 N = parentcount, id =  PTID,
                                 data = booldata,
                                 treatment = treatment,
                                 weights = NULL,
                                 rate = 1, updateLag = 10,
                                 nsamp = 100,
                                 centerCovariance = TRUE,
                                 maxIter = 25, tol = 1e-03)
  data <- booldata
  posteriors <- fit$posteriors[, 3]
  populations <- unique(data$Subset)
  plotList <- lapply(1:length(selected_populations), function(x) x)
  for(i in 1:length(populations)) {
    tempdat <- subset(data, Subset == populations[i])
    prop <- log(tempdat$count / tempdat$parentcount)
    stimprop <- prop[tempdat$stim == "stim"]
    ctrlprop <- prop[tempdat$stim != "stim"]
    vaccine <- (tempdat$vaccine == "VACCINE")[tempdat$stim == "stim"]
    plotList[[i]] <- data.frame(vaccine = vaccine, stimprop = stimprop,
                                ctrlprop = ctrlprop, population = populations[i],
                                posterior = posteriors)
  }
  forPlot <- do.call("rbind", plotList)

  #forPlot <- forPlot[-which(forPlot$posterior > 1), ]
  require(ggplot2)
  try(print(ggplot(forPlot) +
    geom_point(aes(x = ctrlprop, y = stimprop, col = posterior, shape = !vaccine), fill = "white") +
    theme_bw() + geom_abline(intercept = 0, slope = 1) +
    scale_colour_gradientn(colours=rainbow(4)) +
    facet_wrap(~ population, ncol = 5)))

  par(mfrow = c(1, 2))
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
       #main = populations[i])
  abline(a = 0, b = 1)
  lines(uniqueNominal, power, col = "blue", lty = 2)
  legend("topright", col = c("red", "blue"), lty = 1:2, legend = c("FDR", "Power"))
  abline(v = c(.01, .05, .1), h = c(.75, .9, .95), col = "grey")
  rocfit <- roc(vaccine ~ posteriors)
  print(plot(rocfit, main = round(rocfit$auc, 3)))
}
dev.off()
