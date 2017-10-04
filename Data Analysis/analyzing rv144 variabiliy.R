# Preparing data -----------------
require(pROC)
require(reshape2)
load("data/rv144_booleans.rda")
bySubset <- by(data.frame(booleans$stim, booleans$nonstim), booleans$Subset, function(x) x)
largerThanThershold <- sapply(bySubset, function(x) colSums(x >5))

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
groups <- list(c(1:38), c(1:20), c(21:38))
resList <- list()
for(i in 1:length(groups)) {
  temppost <- postList[groups[[i]]]
  aucs <- matrix(nrow = length(temppost), ncol = ncol(fit$posteriors) - 1)
  for(j in 1:length(temppost)) {
    aucs[j, ] <- apply(temppost[[j]], 2, function(x) roc(infect ~ x)$auc)
  }
  means <- colMeans(aucs)
  sds <- apply(aucs, 2, sd) /sqrt(nrow(aucs))
  resList[[i]] <- data.frame(lci = means - 2 * sds, estimate = means, uci = means + 2 * sds)
  rownames(resList[[i]]) <- colnames(fit$posteriors)[-1]
  resList[[i]] <- resList[[i]][order(resList[[i]][, 2], decreasing = TRUE), ]
}


# Analysis for a single fit ----------------------
reps <- 200
aucs <- matrix(nrow = reps, ncol = ncol(fit$posteriors) - 1)
rprob <- matrix(nrow = reps, ncol = ncol(fit$posteriors) - 1)
niters <- nrow(subjList[[1]])
for(i in 1:reps) {
  samp <- sample.int(niters, niters, replace = TRUE)
  post <- t(sapply(subjList, function(x) colMeans(x[samp, ])))
  aucs[i, ] <- apply(post, 2, function(x) roc(infect ~ x)$auc)
  rprob[i, ] <- colMeans(post)
  cat(i, " ")
}

centerAUC <- apply(fit$posteriors[, -1], 2, function(x) roc(infect ~ x)$auc)
centerRprob <- colMeans(fit$posteriors[, -1])
CI <- matrix(ncol = 2, nrow = length(centerAUC))
probCI <- matrix(ncol = 2, nrow = length(centerAUC))
for(i in 1:length(centerAUC)) {
  CI[i, ] <- quantile(centerAUC[i] - (aucs[, i] - mean(aucs[, i])), probs = c(0.025, 0.975))
  probCI[i, ] <- quantile(centerRprob[i] - (rprob[, i] - mean(rprob[, i])), probs = c(0.025, 0.975))
}
CI <- cbind(CI[, 1], centerAUC, CI[, 2])
probCI <- cbind(probCI[, 1], centerRprob, probCI[, 2])
CI <- data.frame(CI)
probCI <- data.frame(probCI)
names(CI) <- c("lowerCI", "estimate", "upperCI")
names(probCI) <- c("lowerCI", "estimate", "upperCI")

colMeans(aucs)
max(colMeans(aucs))
sds <- apply(aucs, 2, sd)
names(sds) <- names(fit$posteriors)[-1]
