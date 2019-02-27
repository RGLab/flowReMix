library(flowReMix)
library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)

# Parameters --------
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)
ncores <- 1
mode <- "M0M3children"
mode <- "M3children"
mode <- "M0M3stimTreatment"
mode <- "M3infantsChildren"
mode <- "M0M3stimInfants"
mode <- "M0M3all"
mode <- "M0M3allB"

# Loading Data ---------
tempdat <- readRDS(file = "data/MAL067_screened_C.rds")
tempdat$stim <- factor(tempdat$stim)
if(mode == "M0M3children") {
  tempdat$treatment <- factor(with(tempdat, visitno == "M3" & stim != "ctrl"))
  tempdat <- subset(tempdat, agec == "5-17m")
  formula <- formula(cbind(count, parentcount - count) ~ stim  + treatment)
} else if(mode == "M0M3stimTreatment") {
  tempdat$treatment <- tempdat$stim
  tempdat$interaction <- factor(with(tempdat, visitno == "M3" & stim != "ctrl"))
  tempdat <- subset(tempdat, agec == "5-17m")
  formula <- formula(cbind(count, parentcount - count) ~ treatment * visitno)
} else if(mode == "M3children") {
  tempdat$treatment <- tempdat$stim
  tempdat <- subset(tempdat, agec == "5-17m")
  tempdat <- subset(tempdat, visitno == "M3")
  formula(cbind(count, parentcount - count) ~ treatment)
} else if(mode == "M0M3all") {
  tempdat$treatment <- factor(with(tempdat, visitno == "M3" & stim != "ctrl"))
  formula <- formula(cbind(count, parentcount - count) ~ treatment + stim + agec)
} else if(mode == "M0M3allB") {
  tempdat$treatment <- factor(with(tempdat, visitno == "M3" & stim != "ctrl"))
  formula <- formula(cbind(count, parentcount - count) ~ treatment + stim * sagec)
} else if(mode == "M0M3stimInfants") {
  tempdat$treatment <- tempdat$stim
  formula <- formula(cbind(count, parentcount - count) ~ treatment*agec*visitno)
} else if(mode == "M3infantsChildren"){
  tempdat$treatment <- tempdat$stim
  tempdat <- subset(tempdat, visitno == "M3")
  formula <- formula(cbind(count, parentcount - count) ~ treatment*agec)
}

# Analysis Parameters -----
configurations <- expand.grid(iterations = c(60),
                              mcEM = TRUE,
                              disp = c(50),
                              npost = c(1),
                              seed = 101:200)
config <- configurations[setting, ]
niter <- config[["iterations"]]
lag <- round(niter / 3)
mcEM <- config[["mcEM"]]
npost <- config[["npost"]]
seed <- config[["seed"]]
disp <- config[["disp"]]
stabilityReps <- 50

# Control Object --------
control <- flowReMix_control(updateLag = lag, nsamp = 50, initMHcoef = 1,
                             keepEach = 5,
                             nPosteriors = npost, centerCovariance = TRUE,
                             maxDispersion = disp * 1000, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             lastSample = NULL, isingInit = -log(99),
                             markovChainEM = mcEM,
                             initMethod = "robust",
                             preAssignCoefs = 1,
                             seed = seed,
                             ncores = ncores,
                             isingWprior = FALSE,
                             zeroPosteriorProbs = FALSE,
                             isingStabilityReps = stabilityReps,
                             randStabilityReps = 0,
                             learningRate = 0.75,
                             keepWeightPercent = 0.9,
                             sampleNew = FALSE,
                             stabilityAND = FALSE, stabilityGamma = 0)

# Analysis --------
tempdat$visitno <- factor(tempdat$visitno, levels = c("M0", "M3"))
tempdat$agec <- factor(tempdat$agec)
tempdat$subset <- factor(tempdat$subset)
# tempdat$interaction <- tempdat$visitno:tempdat$stim
fit <- flowReMix(cbind(count, parentcount - count) ~ treatment + stim * agec,
                 subject_id = ptid,
                 cell_type = subset,
                 cluster_variable = treatment,
                 data = tempdat,
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "robust",
                 cluster_assignment = TRUE,
                 iterations = niter,
                 parallel = TRUE,
                 verbose = TRUE, control = control)

# file <- paste("results/malVaccine_newInterM3infantsStim_A_",
#               "disp", disp,
#               "seed", seed,
#               "npost", npost,
#               "niter", niter,
#               ".rds", sep ="")
file <- paste("results/malVaccine_newInterM0M3infantsB_B_",
              "disp", disp,
              "seed", seed,
              "npost", npost,
              "niter", niter,
              ".rds", sep ="")
saveRDS(fit, file = file)


