cpus <- 1
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
                              seed = c(11:40),
                              prior = c(0),
                              wprior = c(FALSE, TRUE),
                              method = c("SA", "MC"))
configurations$method <- as.character(configurations$method)
config <- configurations[setting, ]
print(config)
niter <- (config[["niters"]])
seed <- (config[["seed"]])
prior <- (config[["prior"]])
method <- config["method"][[1]]
wprior <- (config[["wprior"]])
if(method == "SA") {
  npost <- 1
  lag <- 10
  mcEM <- FALSE
} else {
  npost <- 3
  lag <- 20
  mcEM <- TRUE
}


# Analysis -------------
library(flowReMix)
control <- flowReMix_control(updateLag = lag, nsamp = 200,
                             keepEach = 10, initMHcoef = 2.5,
                             nPosteriors = npost, centerCovariance = FALSE,
                             maxDispersion = 10^4, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             isingInit = -log(99),
                             seed = seed,
                             ncores = cpus, preAssignCoefs = 1,
                             prior = prior, isingWprior = wprior,
                             markovChainEM = mcEM,
                             initMethod = "robust",
                             learningRate = 0.6, keepWeightPercent = .9,
                             isingStabilityReps = 100)

booldata$subset <- factor(booldata$subset)
system.time(fit <- flowReMix(cbind(count, parentcount - count) ~ treatment,
                             subject_id = ptid,
                             cell_type = subset,
                             cluster_variable = treatment,
                             data = booldata,
                             covariance = "sparse",
                             ising_model = "sparse",
                             regression_method = "robust",
                             iterations =  niter,
                             cluster_assignment = TRUE,
                             parallel = TRUE, keepSamples = FALSE,
                             verbose = TRUE, control = control))
file <- paste("results/rv144_32_niter", niter, "npost", npost, "seed", seed, "prior", prior, "pising",as.numeric(wprior), method,".rds", sep = "")
print(file)
saveRDS(fit, file = file)
