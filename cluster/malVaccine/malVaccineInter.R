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

# Loading Data ---------
tempdat <- readRDS(file = "data/MAL067_screened_B.rds")
tempdat$stim <- factor(tempdat$stim)
tempdat$interaction <- with(tempdat, visitno == "M3" & stim != "ctrl") %>% factor()
# obsPerSubject <- tempdat %>% select(ptid, visitno, stim) %>% unique() %>%
#   group_by(ptid) %>%
#   summarize(nObs = length(ptid))
# table(obsPerSubject$nObs)
# obsPerSubject <- subset(obsPerSubject, nObs == 4)
# tempdat <- subset(tempdat, ptid %in% obsPerSubject$ptid)

# Analysis Parameters -----
configurations <- expand.grid(iterations = c(60),
                              mcEM = TRUE,
                              disp = c(50),
                              npost = c(1),
                              seed = 1:50)
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
# tempdat <- subset(tempdat, !(population %in% c("CD154+", "GzB+")))
tempdat$subset <- factor(tempdat$subset)
tempdat$isInfant <- tempdat$agec == "6-12w"
tempdat$isChild <- tempdat$agec == "5-17m"
fit <- flowReMix(cbind(count, parentcount - count) ~ agec*stim + stim*visitno,
                 subject_id = ptid,
                 cell_type = stimCellType,
                 cluster_variable = stim,
                 data = tempdat,
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "robust",
                 cluster_assignment = TRUE,
                 iterations = niter,
                 parallel = TRUE,
                 verbose = TRUE, control = control)

file <- paste("results/malVaccine_stimInfants_F_",
              "disp", disp,
              "seed", seed,
              "npost", npost,
              "niter", niter,
              ".rds", sep ="")
saveRDS(fit, file = file)

# saveRDS(fit, "data analysis/results/malVaccine_tempfit.rds")
fit$data$malaria[fit$data$vaccine2c == "comparator"] <- NA
summary(fit)
