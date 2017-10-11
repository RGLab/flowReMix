# Processing dataset -------------------
getExpression <- function(str) {
  first <- substr(str, 1, 7)
  second <- substr(str, 8, nchar(str))
  second <- strsplit(second, "")[[1]]
  seperators <- c(0, which(second %in% c("-", "+")))
  expressed <- list()
  for(i in 2:length(seperators)) {
    if(second[seperators[i]] == "+") {
      expressed[[i]] <- paste(second[(seperators[(i - 1)] + 1) : seperators[i]], collapse = '')
    }
  }

  expressed <- paste(unlist(expressed), collapse = '')
  expressed <- paste(first, expressed, sep = '')
  return(expressed)
}

assign <- function(x) {
  x$prop <- x$count / x$parentcount
  assign <- as.numeric(by(x, x$subset, function(y) max(y$prop[y$stim != 0]) > min(y$prop[y$stim == 0])))
  assign[assign == 1] <- -1
  result <- data.frame(ptid = x$ptid[1], subset = unique(x$subset), assign = assign)
  return(result)
}

# Loading Data --------------------------------
# hvtn <- read.csv(file = "data/merged_505_stats.csv")
# names(hvtn) <- tolower(names(hvtn))
# hvtn <- subset(hvtn, !is.na(ptid))
# saveRDS(hvtn, file = "data/505_stats.rds")


# Getting Demographic data ------------------------
demo <- read.csv(file = "data/primary505.csv")
infect <- data.frame(ptid = demo$ptid, status = demo$HIVwk28preunbl)
infect <- subset(infect, infect$ptid %in% hvtn$ptid)

# Getting marginals -----------------------------
library(flowReMix)
hvtn <- readRDS(file = "data/505_stats.rds")
length(unique(hvtn$name))
length(unique(hvtn$ptid))
length(unique(hvtn$population))
unique(hvtn$population)
unique(hvtn$stim)
nchars <- nchar(as.character(unique(hvtn$population)))
#marginals <- unique(hvtn$population)[nchars < 26]
marginals <- unique(hvtn$population)[nchars == 26]
marginals <- subset(hvtn, population %in% marginals)
marginals <- subset(marginals, stim %in% c("negctrl", "VRC ENV A",
                                           "VRC ENV B", "VRC ENV C",
                                           "VRC GAG B", "VRC NEF B",
                                           "VRC POL 1 B", "VRC POL 2 B"))
marginals <- subset(marginals, !(population %in% c("4+", "8+")))
marginals <- subset(marginals, !(population %in% c("8+/107a-154-IFNg-IL2-TNFa-", "4+/107a-154-IFNg-IL2-TNFa-")))
marginals$stim <- factor(as.character(marginals$stim))
marginals$population <- factor(as.character(marginals$population))

# Descriptives -------------------------------------
library(ggplot2)
marginals$prop <- marginals$count / marginals$parentcount
# ggplot(marginals) + geom_boxplot(aes(x = population, y = log(prop), col = stim))

require(dplyr)
negctrl <- subset(marginals, stim == "negctrl")
negctrl <- summarize(group_by(negctrl, ptid, population), negprop = mean(prop))
negctrl <- as.data.frame(negctrl)
marginals <- merge(marginals, negctrl, all.x = TRUE)

# ggplot(subset(marginals, stim != "negctrl" & parent == "4+")) +
#   geom_point(aes(x = log(negprop), y = log(prop)), size = 0.25) +
#   facet_grid(stim ~ population, scales = "free") +
#   theme_bw() +
#   geom_abline(intercept = 0, slope = 1)

# Setting up data for analysis ---------------------------
unique(marginals$stim)
gag <- subset(marginals, stim %in% c("VRC GAG B", "negctrl"))
gag$subset <- factor(paste("gag", gag$population, sep = "/"))
gag$stimGroup <- "gag"
pol <-subset(marginals, stim %in% c("negctrl", "VRC POL 1 B", "VRC POL 2 B"))
pol$subset <- factor(paste("pol", pol$population, sep = "/"))
pol$stimGroup <- "pol"
env <- subset(marginals, stim %in% c("negctrl", "VRC ENV C", "VRC ENV B", "VRC ENV A"))
env$subset <- factor(paste("env", env$population, sep = "/"))
env$stimGroup <- "env"
nef <- subset(marginals, stim %in% c("negctrl", "VRC NEF B"))
nef$subset <- factor(paste("nef", nef$population, sep = "/"))
nef$stimGroup <- "nef"
subsetDat <- rbind(gag, pol, env, nef)
subsetDat$stim <- as.character(subsetDat$stim)
subsetDat$stim[subsetDat$stim == "negctrl"] <- 0
subsetDat$stim <- factor(subsetDat$stim)

# Converting subset names ------------------
subsets <- as.character(unique(subsetDat$subset))
expressed <- sapply(subsets, getExpression)
map <- cbind(subsets, expressed)
subsetDat$subset <- as.character(subsetDat$subset)
for(i in 1:nrow(map)) {
  subsetDat$subset[which(subsetDat$subset == map[i, 1])] <- map[i, 2]
}
subsetDat$subset <- factor(subsetDat$subset)

# Getting outcomes -------------------------------
treatmentdat <- read.csv(file = "data/rx_v2.csv")
names(treatmentdat) <- tolower(names(treatmentdat))
treatmentdat$ptid <- factor(gsub("-", "", (treatmentdat$ptid)))
treatmentdat <- subset(treatmentdat, ptid %in% unique(subsetDat$ptid))

# Finding problematic subsets?
keep <- by(subsetDat, list(subsetDat$subset), function(x) mean(x$count > 1) > 0.02)
keep <- names(keep[sapply(keep, function(x) x)])
#result$subsets[result$qvals < 0.1] %in% keep
subsetDat <- subset(subsetDat, subset %in% keep)
subsetDat$subset <- factor(as.character(subsetDat$subset))
subsetDat <- merge(subsetDat, treatmentdat)
subsetDat <- merge(subsetDat, infect)
subsetDat$hiv <- subsetDat$status
subsetDat$hiv[subsetDat$control == 1] <- NA
subsetDat$vaccine <- 1 - subsetDat$control
subsetDat$batch <- factor(subsetDat$batch..)
subsetDat$stimGroup <- factor(subsetDat$stimGroup)

# Loading files ----------------------

filenames <- as.list(dir(path = 'data analysis/results', pattern="hvtn_2_*"))[1:88]
filenames <- lapply(filenames, function(x) paste0('data analysis/results/', x))[-c(3, 4)]
post <- list()
postList <- list()
rm(temp)
for(i in 1:length(filenames)) {
  load(file = filenames[[i]])
  try(fit <- temp)
  post[[i]] <- fit$posteriors[, -1]
  postList[[i]] <- fit$posteriors[, -1]
}
post <- Reduce("+", post) / length(filenames)
# fit$data <- booldata
fit$posteriors[, -1] <- post

# Setting up target variables ------
hiv <- data.frame(ptid = subsetDat$ptid, hiv = subsetDat$hiv, vaccine = subsetDat$vaccine)
hiv <- unique(hiv)
postid <- fit$posteriors[, 1, drop = FALSE]
hiv <- merge(postid, hiv)

# Variability at the fit level --------------------
groups <- list(c(1:47), c(48:88), c(1:88))
groups <- list(c(1:47), c(48:88))
resList <- list()
rpList <- list()
auclist <- list()
problist <- list()
library(pROC)
for(i in 1:length(groups)) {
  temppost <- postList[groups[[i]]]
  aucs <- matrix(nrow = length(temppost), ncol = ncol(fit$posteriors) - 1)
  responseProb <- matrix(nrow = length(temppost), ncol = ncol(fit$posteriors) - 1)
  for(j in 1:length(temppost)) {
    try(aucs[j, ] <- apply(temppost[[j]], 2, function(x) roc(hiv[, 2] ~ x)$auc))
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

  colnames(aucs) <- colnames(fit$posteriors)[-1]
  aucs <- data.frame(aucs)
  auclist[[i]] <- aucs
  colnames(responseProb) <- colnames(aucs)
  problist[[i]] <- data.frame(responseProb)
  print(cbind(resList[[i]], rpList[[i]]))
}

library(dplyr)
library(reshape2)
aucplot <- auclist
aucplot[[1]]$iterations <- 30
aucplot[[2]]$iterations <- 60
aucplot <- do.call("rbind", aucplot)
aucplot <- melt(aucplot, id = "iterations")
names(aucplot)[2:3] <- c("subset", "auc")
ggplot(aucplot) + geom_boxplot(aes(x = subset, y = auc, col = factor(iterations))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

probplot <- problist
probplot[[1]]$iterations <- 30
probplot[[2]]$iterations <- 60
probplot <- do.call("rbind", probplot)
probplot <- melt(probplot, id = "iterations")
names(probplot)[2:3] <- c("subset", "responseProb")
ggplot(probplot) +
  geom_boxplot(aes(x = subset, y = responseProb, col = factor(iterations), outlier.size = 0.01)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Graphical models -----------------
filenames <- as.list(dir(path = 'data analysis/results', pattern="stab_hvtn_2_*"))[1:47]
filenames <- as.list(dir(path = 'data analysis/results', pattern="stab_hvtn_2_*"))[48:88]
filenames <- lapply(filenames, function(x) paste0('data analysis/results/', x))[-c(3, 4)]
net <- list()
for(i in 1:length(filenames)) {
  load(file = filenames[[i]])
  net[[i]] <- stability$network
  # net[[i]] <- stab$network
}
agg <- matrix(0, ncol = ncol(net[[1]]), nrow = nrow(net[[1]]))
rownames(agg) <- rownames(net[[1]])
colnames(agg) <- colnames(net[[1]])
for(i in 1:nrow(agg)) {
  for(j in 1:ncol(agg)) {
    agg[i, j] <- median(sapply(net, function(x) x[i, j]))
  }
}
net <- Reduce("+", net) #/ length(net)
stability$network <- agg
stab$network <- net
plot(stability, threshold = .8, fill = rocResults$auc)
