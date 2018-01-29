cpus <- 4
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)

# Load data ------
booldata <- readRDS(file = "data/rv144bool.rds")

# Configurations --------------------
configurations <- expand.grid(niters = c(60),
                              npost = c(5),
                              seed = c(1:100))
config <- configurations[setting, ]
niter <- config[[1]]
npost <- config[[2]]
seed <- config[[3]]

# Analysis -------------
library(flowReMix)
control <- flowReMix_control(updateLag = round(niter / 3), nsamp = 50,
                             keepEach = 5, initMHcoef = 2.5,
                             nPosteriors = npost, centerCovariance = FALSE,
                             maxDispersion = 10^4, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             isingInit = -log(99),
                             seed = seed,
                             ncores = cpus, preAssignCoefs = 1,
                             prior = 0, isingWprior = FALSE,
                             markovChainEM = TRUE,
                             initMethod = "robust",
                             learningRate = 0.6, keepWeightPercent = .9,
                             lastSample = NULL,
                             isingStabilityReps = 0)

sampid <- function(data, vaccine == TRUE) {
  a <- subset(data, sample(data$ptid[data$], 1) == ptid)
  a$ptid <- runif(1)
  return(a)
}
booldata$subset <- factor(booldata$subset)
bootsamp <- do.call("rbind", replicate(length(unique(booldata$ptid)), sampid(booldata), simplify = FALSE))
bootsamp$ptid <- factor(bootsamp$ptid)
system.time(fit <- flowReMix(cbind(count, parentcount - count) ~ treatment,
                             subject_id = ptid,
                             cell_type = subset,
                             cluster_variable = treatment,
                             data = booldata,
                             covariance = "diagonal",
                             ising_model = "none",
                             regression_method = "robust",
                             iterations =  niter,
                             cluster_assignment = preAssignment,
                             parallel = TRUE, keepSamples = FALSE,
                             verbose = TRUE, control = control))
file <- paste("results/rv144_42_boot", niter, "npost", npost, "seed", seed, "MC.rds", sep = "")
saveRDS(fit, file = file)
