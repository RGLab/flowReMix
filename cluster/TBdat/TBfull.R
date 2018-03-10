library(magrittr)
library(flowReMix)
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)
ncores <- 4

assign <- function(x) {
  x$prop <- x$count / x$parentcount
  assign <- as.numeric(by(x, x$subset, function(y) max(y$prop[y$stim != "aUNS"]) > min(y$prop[y$stim == "aUNS"])))
  assign[assign == 1] <- -1
  result <- data.frame(ptid = x$ptid[1], subset = unique(x$subset), assign = assign)
  return(result)
}

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

# Defining subsets wo stim --------------
tbdat$subset <- interaction(tbdat$parent, tbdat$population, sep = "/")
tbdat <- subset(tbdat, stim != "EBV")
tbdat$stim[tbdat$stim == "UNS"] <- "ctrl"
tbdat$stim <- factor(tbdat$stim, levels = c("ctrl", "MP", "Mtbaux", "P1", "P2", "P3"))
jointcounts <- by(tbdat, tbdat$subset, function(x) mean(x$count > 0))
jointcounts <- data.frame(names(jointcounts), as.numeric(unlist(jointcounts)))
tbdat <- data.frame(tbdat)

# Defining subsets w stim ---------------------
stimDat <- stimulationModel(tbdat, subset, stim,
                            controls = "ctrl",
                            stim_groups = list(MP = c("MP", "P1", "P2", "P3"),
                                               Mtbaux = "Mtbaux"))
stimDat$subset <- stimDat$stimCellType
stimDat$stimCellType <- NULL
stimcounts <- by(stimDat, stimDat$subset, function(x) mean(x$count > 0))
stimcounts <- data.frame(names(stimcounts), as.numeric(unlist(stimcounts)))

# Choosing subset of data for analysis -----------------
jointkeep <- jointcounts[jointcounts[, 2] >= 0.2, 1]
stimkeep <- stimcounts[stimcounts[, 2] >= 0.2, 1]
tbdat <- subset(tbdat, subset %in% jointkeep)
stimDat <- subset(stimDat, subset %in% stimkeep)
tbdat$subset <- factor(tbdat$subset)
stimDat$subset <- factor(stimDat$subset)

# Analysis Setting -------------
configurations <- expand.grid(mcEM = c(TRUE),
                              maxdisp = c(10, 50, 100),
                              seed = 1:50,
                              npost = c(10),
                              niter = c(60))

mcEM <- configurations[["mcEM"]][setting]
seed <- configurations[["seed"]][setting]
npost <- configurations[["npost"]][setting]
niter <- configurations[["niter"]][setting]
maxdisp <- configurations[["maxdisp"]][setting]
lag <- round(niter / 3)

# Analysis ----------------------------------
control <- flowReMix_control(updateLag = lag, nsamp = 50, initMHcoef = 1,
                             keepEach = 5,
                             nPosteriors = npost, centerCovariance = TRUE,
                             maxDispersion = maxdisp * 1000, minDispersion = 10^8,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             lastSample = NULL, isingInit = -7,
                             markovChainEM = mcEM,
                             initMethod = "robust",
                             preAssignCoefs = 1,
                             seed = seed,
                             ncores = ncores,
                             prior = -0.5,
                             isingWprior = FALSE,
                             zeroPosteriorProbs = FALSE,
                             isingStabilityReps = 100,
                             randStabilityReps = 0,
                             learningRate = 0.75,
                             keepWeightPercent = 0.9,
                             sampleNew = FALSE)

fit <- flowReMix(cbind(count, parentcount - count) ~ stim,
                 subject_id = ptid,
                 cell_type = subset,
                 cluster_variable = stim,
                 data = stimDat,
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "robust",
                 cluster_assignment = TRUE,
                 iterations = niter,
                 parallel = TRUE,
                 verbose = TRUE, control = control)

file <- paste("results/TBdat_rob_A_",
              "disp", maxdisp,
              "seed", seed,
              "npost", npost,
              "niter", niter,
              ".rds", sep ="")
saveRDS(fit, file = file)





