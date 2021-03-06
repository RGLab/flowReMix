cpus <- 2
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)

assign <- function(x) {
  x$prop <- x$count / x$parentcount
  assign <- as.numeric(by(x, x$subset, function(y) y$prop[1] > y$prop[2]))
  assign[assign == 1] <- -1
  result <- data.frame(ptid = x$ptid[1], subset = unique(x$subset), assign = assign)
  return(result)
}

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
names(booldata) <- tolower(names(booldata))

# Configurations --------------------
configurations <- expand.grid(niters = c(40, 80),
                              npost = c(1),
                              seed = c(1:50),
                              method = c("MC", "SA"))
config <- configurations[setting, ]
niter <- config[[1]]
npost <- config[[2]]
seed <- config[[3]]
method <- config[[4]]
if(method == "MC") {
  lag <- round(niter / 2)
  keepEach <- 10
  mcem <- TRUE
} else if(method == "SA") {
  lag <- 10
  keepEach <- 20
  mcem <- FALSE
}

# Analysis -------------
library(flowReMix)
control <- flowReMix_control(updateLag = lag, nsamp = 100,
                             keepEach = keepEach, initMHcoef = 2.5,
                             nPosteriors = npost, centerCovariance = FALSE,
                             maxDispersion = 10^3, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 100,
                             isingInit = -log(999),
                             seed = seed,
                             ncores = cpus, preAssignCoefs = 1,
                             prior = 0, isingWprior = FALSE,
                             markovChainEM = mcem,
                             initMethod = "robust",
                             learningRate = 0.6, keepWeightPercent = .9,
                             isingStabilityReps = 200)

booldata$subset <- factor(booldata$subset)
preAssignment <- do.call("rbind", by(booldata, booldata$ptid, assign))
system.time(fit <- flowReMix(cbind(count, parentcount - count) ~ treatment,
                             subject_id = ptid,
                             cell_type = subset,
                             cluster_variable = treatment,
                             data = booldata,
                             covariance = "sparse",
                             ising_model = "sparse",
                             regression_method = "robust",
                             iterations =  niter,
                             cluster_assignment = preAssignment,
                             parallel = TRUE, keepSamples = FALSE,
                             verbose = TRUE, control = control, newSampler = TRUE))
file <- paste("results/rv144_29_niter", niter, "npost", npost, "seed", seed, method,".rds", sep = "")
saveRDS(fit, file = file)
