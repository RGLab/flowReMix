# Preparing data -----------------
require(pROC)
require(reshape2)
load("data/rv144_booleans.rda")
bySubset <- by(data.frame(booleans$stim, booleans$nonstim), booleans$Subset, function(x) x)
largerThanThershold <- sapply(bySubset, function(x) colSums(x > 5))

booldata <- melt(booleans, c("PTID", "Subset"))
names(booldata)[3:4] <- c("stim", "count")

booldata <- by(booldata, INDICES = list(booldata$PTID, booldata$stim), function(x) {
  x$parentcount <- sum(x$count)
  return(x)
})
booldata <- do.call("rbind", booldata)

booldata <- subset(booldata, Subset != "!TNFa&!IFNg&!IL4&!IL2&!CD154&!IL17a")
booldata$treatment <- as.numeric(booldata$stim == "stim")
uniquepop <- unique(booldata$Subset)
booldata <- with(booldata, booldata[order(Subset, PTID, stim, decreasing = FALSE), ])
booldata <- subset(booldata, !is.na(Subset))
allsubset <- booldata
booldata <- with(booldata, booldata[order(Subset, PTID, stim, decreasing = FALSE), ])

# Naming ------------------
subsets <- unique(booldata$Subset)
booldata$Subset <- as.character(booldata$Subset)
nfunctions <- numeric(length(subsets))
for(i in 1:length(subsets)) {
  split <- strsplit(as.character(subsets[i]), "&")[[1]]
  first <- substr(split, 1, 1)
  nfunction <- sum(first != "!")
  nfunctions[i] <- nfunction
  name <- paste(split[first != "!"], collapse = ",")
  booldata$nfunction[booldata$Subset == subsets[[i]]] <- nfunction
  booldata$Subset[booldata$Subset == subsets[[i]]] <- name
}
subsets <- unique(booldata$Subset)
booldata <- with(booldata, booldata[order(Subset, PTID, stim, decreasing = FALSE), ])
# booldata <- merge(booldata, data.frame(ptid = correlates$ptid,
#                                        IgAprim = correlates$IgAprim,
#                                        V2prim = correlates$V2prim))
names(booldata) <- tolower(names(booldata))


# Getting vaccine information --------------------
data("rv144")
rv144 <- rv144[order(rv144$ptid), ]
vaccine <- (by(rv144, rv144$ptid, function(x) x$vaccine[1] == "VACCINE"))
vaccine <- data.frame(ptid = names(vaccine), vaccine = as.numeric(vaccine))
vaccinemat <- vaccine[vaccine$ptid %in% booldata$ptid, ]

# Getting infection status
data("rv144_correlates_data")
correlates <- rv144_correlates_data
correlates <- correlates[order(as.character(correlates$PTID)), ]
infection <- correlates$infect.y
names(correlates)[1] <- "ptid"
vaccine <- correlates[, c(1, 62, 64)]
vaccine$ptid <- factor(vaccine$ptid, levels = levels(booldata$ptid))
booldata <- with(booldata, booldata[order(subset, ptid, stim, decreasing = FALSE), ])
booldata <- merge(booldata, vaccine, sort = FALSE)
booldata$hiv <- NA
booldata$hiv[booldata$infect.y == "INFECTED"] <- TRUE
booldata$hiv[booldata$infect.y == "NON-INFECTED"] <- FALSE
booldata <- with(booldata, booldata[order(subset, ptid, stim, decreasing = FALSE), ])
booldata$subset <- factor(booldata$subset)
booldata$stim <- factor(booldata$stim, levels = c("nonstim", "stim"))

# Getting result files --------------------------
filenames <- as.list(dir(path = 'data analysis/results', pattern="rv144_20__*"))
filenames <- as.list(dir(path = 'data analysis/results', pattern="rv144_23__*"))[1:197]
filenames <- as.list(dir(path = 'data analysis/results', pattern="rv144_24__*"))[1:197]
filenames <- as.list(dir(path = 'data analysis/results', pattern="rv144_25__*"))
filenames <- lapply(filenames, function(x) paste0('data analysis/results/', x))[-c(3, 4)]
post <- list()
postList <- list()
for(i in 1:length(filenames)) {
  load(file = filenames[[i]])
  post[[i]] <- fit$posteriors[, -1]
  postList[[i]] <- fit$posteriors[, -1]
}
post <- Reduce("+", post) / length(filenames)
fit$data <- booldata
fit$posteriors[, -1] <- post

# Setting up target variables -----------------
ids <- fit$posteriors[, 1:2]
vaccine[, 1] <- as.character(vaccine[, 1])
vaccine[, 1] <- factor(vaccine[, 1], levels = levels(ids[, 1]))
vaccine <- vaccine[!is.na(vaccine[, 1]), ]
vaccine <- vaccine[order(vaccine[, 1]), ]
ids <- merge(ids, vaccine, all.x = TRUE, all.y = FALSE,
             by = "ptid", sort = FALSE)
vaccination <- ids[, 3]
infectDat <- data.frame(ptid = rv144_correlates_data$PTID, infect = rv144_correlates_data$infect.y)
datId <- as.character(fit$posteriors$ptid)
infectID <- as.character(infectDat$ptid)
infectDat <- infectDat[infectID %in% datId, ]
infectDat$ptid <- factor(as.character(infectDat$ptid), levels = levels(booldata$ptid))
infectDat <- infectDat[order(infectDat$ptid), ]
ids <- merge(ids, infectDat, sort = FALSE)
infect <- ids[, 4]
infect[infect == "PLACEBO"] <- NA
infect <- factor(as.character(infect), levels = c("INFECTED", "NON-INFECTED"))


# Bootstrapping -------------------------
# groups <- list(c(1:38), c(1:20), c(21:38))
groups <- list(c(181:280), c(281:378), c(1:93), c(94:180))
# groups <- list(c(1:98))
resList <- list()
rpList <- list()
for(i in 1:length(groups)) {
  temppost <- postList[groups[[i]]]
  aucs <- matrix(nrow = length(temppost), ncol = ncol(fit$posteriors) - 1)
  responseProb <- matrix(nrow = length(temppost), ncol = ncol(fit$posteriors) - 1)
  for(j in 1:length(temppost)) {
    try(aucs[j, ] <- apply(temppost[[j]], 2, function(x) roc(infect ~ x)$auc))
    try(responseProb[j, ] <- colMeans(temppost[[j]]))
  }
  means <- colMeans(aucs)
  sds <- apply(aucs, 2, sd) /sqrt(nrow(aucs))
  quantiles <- t(apply(aucs, 2, quantile, na.rm = TRUE))
  rpQuantiles <- t(apply(responseProb, 2, quantile, na.rm = TRUE))
  rownames(quantiles) <- colnames(fit$posteriors[-1])
  rownames(rpQuantiles) <- colnames(fit$posteriors[-1])
  # resList[[i]] <- data.frame(lci = means - 2 * sds, estimate = means, uci = means + 2 * sds)
  resList[[i]] <- quantiles
  rpList[[i]] <- rpQuantiles
  rpList[[i]] <- rpList[[i]][order(resList[[i]][, 3], decreasing = TRUE), ]
  resList[[i]] <- resList[[i]][order(resList[[i]][, 3], decreasing = TRUE), ]
  print(cbind(resList[[i]], rpList[[i]]))
}

rep1 <- resList[[2]]
rep2 <- resList[[3]]
rep1 <- rep1[order(rep2[, 3], decreasing = TRUE), ]
colnames(rep1) <- paste(50, "iter", colnames(rep1), sep = "")
colnames(rep2) <- paste(100, "iter", colnames(rep2), sep = "")
rep <- cbind(rep1, rep2)
round(rep, 3)

names(rpList) <- paste(c(40, 80, 160, 320), "iterations")

# Bootstrapping the variabilty of an aggregate fit --------------
reps <- 100
results <- list()
presults <- list()
for(g in 1:length(groups)) {
  aucs <- matrix(nrow = reps, ncol = ncol(fit$posteriors) - 1)
  rprobs <- matrix(nrow = reps, ncol = ncol(fit$posteriors) - 1)
  for(r in 1:reps) {
    temppost <- postList[groups[[g]]]
    samp <- temppost[sample.int(length(temppost), length(temppost), replace = TRUE)]
    samp <- Reduce("+", samp) / length(samp)
    aucs[r, ] <- apply(samp, 2, function(x) roc(infect ~ x)$auc)
    rprobs[r, ] <- colMeans(samp)
  }
  means <- colMeans(rprobs)
  sds <- apply(rprobs, 2, sd)
  pres <- data.frame(lci = means - 2 * sds, estimate = means, uci = means + 2 * sds)
  rownames(pres) <- colnames(fit$posteriors)[-1]

  means <- colMeans(aucs)
  sds <- apply(aucs, 2, sd)
  res <- data.frame(lci = means - 2 * sds, estimate = means, uci = means + 2 * sds)
  rownames(res) <- colnames(fit$posteriors)[-1]
  pres <- pres[order(res[, 2], decreasing = TRUE), ]
  res <- res[order(res[, 2], decreasing = TRUE), ]
  print(cbind(res, pres))
  results[[g]] <- res
  presults[[g]] <- pres
}

# temp <- results
# results[[3]] <- results[[3]][order(results[[2]][, 2]), ]
report <- cbind(results[[1]], results[[2]])
names(report) <- c("lci50", "median50", "uci50", "lci100", "median100", "uci100")

# Concensus graph -----------------------
filenames <- as.list(dir(path = 'data analysis/results', pattern="rv144_24__*"))[200:390]
filenames <- as.list(dir(path = 'data analysis/results', pattern="rv144_stab_23__*"))[1:196]
filenames <- as.list(dir(path = 'data analysis/results', pattern="rv144_stab_25__*"))[c(groups[[4]])]
filenames <- lapply(filenames, function(x) paste0('data analysis/results/', x))[-c(3, 4)]
net <- list()
for(i in 1:length(filenames)) {
  load(file = filenames[[i]])
  # net[[i]] <- stability$network
  net[[i]] <- stab$network
}
net <- Reduce("+", net) / length(net)
stability$network <- net
stab$network <- net
plot(stab, threshold = .95, fill = rocResults$auc)

scatter <- plot(fit, type = "scatter", target = vaccine)

# Analysis for a single fit ----------------------
# reps <- 200
# aucs <- matrix(nrow = reps, ncol = ncol(fit$posteriors) - 1)
# rprob <- matrix(nrow = reps, ncol = ncol(fit$posteriors) - 1)
# niters <- nrow(subjList[[1]])
# for(i in 1:reps) {
#   samp <- sample.int(niters, niters, replace = TRUE)
#   post <- t(sapply(subjList, function(x) colMeans(x[samp, ])))
#   aucs[i, ] <- apply(post, 2, function(x) roc(infect ~ x)$auc)
#   rprob[i, ] <- colMeans(post)
#   cat(i, " ")
# }
#
# centerAUC <- apply(fit$posteriors[, -1], 2, function(x) roc(infect ~ x)$auc)
# centerRprob <- colMeans(fit$posteriors[, -1])
# CI <- matrix(ncol = 2, nrow = length(centerAUC))
# probCI <- matrix(ncol = 2, nrow = length(centerAUC))
# for(i in 1:length(centerAUC)) {
#   CI[i, ] <- quantile(centerAUC[i] - (aucs[, i] - mean(aucs[, i])), probs = c(0.025, 0.975))
#   probCI[i, ] <- quantile(centerRprob[i] - (rprob[, i] - mean(rprob[, i])), probs = c(0.025, 0.975))
# }
# CI <- cbind(CI[, 1], centerAUC, CI[, 2])
# probCI <- cbind(probCI[, 1], centerRprob, probCI[, 2])
# CI <- data.frame(CI)
# probCI <- data.frame(probCI)
# names(CI) <- c("lowerCI", "estimate", "upperCI")
# names(probCI) <- c("lowerCI", "estimate", "upperCI")
#
# colMeans(aucs)
# max(colMeans(aucs))
# sds <- apply(aucs, 2, sd)
# names(sds) <- names(fit$posteriors)[-1]
#
