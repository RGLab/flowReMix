ncores <- 15
getwd()

geomean <- function(x) exp(mean(log(x)))
library(flowReMix)
library(pryr)
# Loading data -------------------------
tbdat <- readRDS("data/TB_rozot_booleans.rds")
names(tbdat) <- tolower(names(tbdat))
tbdat$ptid <- sapply(strsplit(tbdat$samplename, "_"), function(x) x[[1]])
tbdat$count[is.na(tbdat$count)] <- 0

# Parsing Expressions ----------------
expressions <- unique(tbdat$pop)
newexp <- character(length(expressions))
for(i in 1:length(expressions)) {
  e <- expressions[i]
  e <- strsplit(e, "&")[[1]]
  res <- c()
  for(j in 1:length(e)) {
    sube <- strsplit(e[j], "-")[[1]][1]
    if(substr(sube, 1, 1) != "!"){
      res <- c(res, sube)
    }
  }
  newexp[i] <- paste(paste(res, collapse = "+"), "+", sep = "")
}
map <- cbind(expressions, newexp)
for(i in 1:nrow(map)) {
  tbdat$population[tbdat$population == map[i, 1]] <- map[i, 2]
}

# Defining subsets w stim ---------------------
tbdat$subset <- paste(tbdat$parent, tbdat$population, sep = "/")
# tbdat <- subset(tbdat, stim != "EBV")
tbdat$stimtemp <- tbdat$stim
# tbdat$stim[tbdat$stim %in% c("P1", "P2", "P3")] <- "MP"
datlist <- list()
stims <- unique(tbdat$stim)
stims <- stims[stims != "UNS"]
for(i in 1:length(stims)) {
  temp <- subset(tbdat, stim %in% c(stims[i], "UNS"))
  temp$stim[temp$stim == "UNS"] <- "ctrl"
  temp$stimgroup <- stims[i]
  datlist[[i]] <- temp
}
tbdat <- do.call("rbind", datlist)
tbdat$stimgroup[tbdat$stim %in% c("P1", "P2", "P3")] <- "MP"
tbdat$subset <- paste(tbdat$stimgroup, tbdat$subset, sep = "/")
nonzerocounts <- by(tbdat, tbdat$subset, function(x) mean(x$count > 0))
nonzerocounts <- data.frame(names(nonzerocounts), as.numeric(unlist(nonzerocounts)))

# Defining subsets wo stim --------------
# tbdat$subset <- paste(tbdat$parent, tbdat$population, sep = "/")
# tbdat <- subset(tbdat, stim != "EBV")
# tbdat$stim[tbdat$stim == "UNS"] <- "ctrl"
# tbdat$stim <- factor(tbdat$stim, levels = c("ctrl", "MP", "Mtbaux", "P1", "P2", "P3"))
# nonzerocounts <- by(tbdat, tbdat$subset, function(x) mean(x$count > 0))
# nonzerocounts <- data.frame(names(nonzerocounts), as.numeric(unlist(nonzerocounts)))

# Choosing subset of data for analysis -----------------
countkeep <- nonzerocounts[nonzerocounts[, 2] >= 0.2, 1]
popkeep <- c("cd4", "MAIT", "NKrainbow", "DCs", "Bcells", "conCD8")
stimkeep <- c("MP", "Mtbaux")
tempdat <- subset(tbdat, (subset %in% countkeep) & (parent %in% popkeep) & (stimgroup %in% stimkeep))
# tempdat <- subset(tbdat, (subset %in% countkeep) & (parent %in% popkeep))
stimlevels <- as.character(unique(tempdat$stim))
stimlevels <- c("ctrl", stimlevels[stimlevels != "ctrl"])
tempdat$stim <- factor(tempdat$stim, levels = stimlevels)
tempdat$subset <- factor(tempdat$subset)
rm(tbdat)

# Visuaizations ------------------------------
library(ggplot2)
tempdat$prop <- tempdat$count / tempdat$parentcount #+ 1 / tempdat$parentcount
subctrl <- subset(tempdat, stim == "ctrl")
substim <- subset(tempdat, stim != "ctrl")
wide <- merge(subctrl, substim, by = c("ptid", "subset"))

# Analysis ----------------------------------
npost <- 30
control <- flowReMix_control(updateLag = 30, nsamp = 50, initMHcoef = 1,
                             nPosteriors = npost, centerCovariance = TRUE,
                             maxDispersion = 10^3, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             ncores = ncores,
                             lastSample = round(100 / npost), isingInit = -log(99),
                             initMethod = "robust")

fit <- flowReMix(cbind(count, parentcount - count) ~ stim,
                 subject_id = ptid,
                 cell_type = subset,
                 cluster_variable = stim,
                 data = tempdat,
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "robust",
                 iterations = 40,
                 parallel = TRUE,
                 verbose = TRUE, control = control)
save(fit, file = "TBsep18iter40post30moreInit.Robj")

stability <- stabilityGraph(fit, type = "ising", cv = FALSE,
                            gamma = 0.25,
                            reps = 200, cpus = ncores)
save(stability, file = "TBising18.Robj")

random <- stabilityGraph(fit, type = "randomEffects", cv = FALSE,
                         gamma = 0.25,
                         reps = 200, cpus = ncores)
save(random, file = "TBrandom18.Robj")






