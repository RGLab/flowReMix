args = as.numeric(commandArgs(trailingOnly = TRUE)[[2]])
settings = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
ncores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))


geomean <- function(x) exp(mean(log(x)))
library(flowReMix)
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
tbdat <- subset(tbdat, stim != "EBV")
tbdat$stimtemp <- tbdat$stim
tbdat$stim[tbdat$stim %in% c("P1", "P2", "P3")] <- "MP"
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
popkeep <- c("cd4", "MAIT", "NKrainbow", "DCs", "Bcells", "conCD*")
stimkeep <- c("P", "MP", "Mtbaux")
# tempdat <- subset(tbdat, (subset %in% countkeep) & (parent %in% popkeep) & (stimgroup %in% stimkeep))
tempdat <- subset(tbdat, (subset %in% countkeep) & (parent %in% popkeep))
#tempdat$stim <- factor(tempdat$stim, levels = c("ctrl", stimkeep))
tempdat$subset <- factor(tempdat$subset)
rm(tbdat)

# Analysis Setting -------------
npost <- 10
niter <- 20
groups <- combn(unique(interaction(tempdat$stimgroup, tempdat$parent)), 2)
tempdat$groups <- interaction(tempdat$stimgroup, tempdat$parent)
instance <- groups[, setting]
tempdat <- subset(tempdat, groups %in% instance)
set.seed(100)
index = (args-1)*1 + settings
seed = round(runif(50, min=1, max=500))

# Analysis ----------------------------------
control <- flowReMix_control(updateLag = 10, nsamp = 50, initMHcoef = 1,
                             keepEach = 5,
                             nPosteriors = npost, centerCovariance = TRUE,
                             maxDispersion = 10^3, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             lastSample = NULL, isingInit = -log(99),
                             markovChainEM = TRUE,
                             initMethod = "robust",
                             preAssignCoefs = 1,
                             seed = as.numeric(seed[index]),
                             ncores = ncores,
                             isingWprior = FALSE,
                             zeroPosteriorProbs = FALSE,
                             isingStabilityReps = 0,
                             randStabilityReps = 0,
                             learningRate = 1,
                             keepWeightPercent = 0.9,
                             sampleNew = FALSE)
tempdat$stim <- tempdat$stimtemp
tempdat$stim[tempdat$stim == "UNS"] <- "aUNS"
tempdat$stim <- as.character(tempdat$stim)
tempdat$stim <- factor(tempdat$stim, levels = sort(unique(tempdat$stim)))
tempdat$subset <- factor(as.character(tempdat$subset))
tempdat <- data.frame(tempdat)
fit <- flowReMix(cbind(count, parentcount - count) ~ stim,
                 subject_id = ptid,
                 cell_type = subset,
                 cluster_variable = stim,
                 data = tempdat,
                 covariance = "diagonal",
                 ising_model = "sparse",
                 regression_method = "robust",
                 cluster_assignment = TRUE,
                 iterations = niter,
                 parallel = TRUE,
                 verbose = FALSE, control = control)

saveRDS(fit, file = paste0("TBfit_",index%%1+1,"_seed_",seed[index],".rds"))
system("sync")
