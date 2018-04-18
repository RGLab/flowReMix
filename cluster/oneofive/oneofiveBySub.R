library(flowReMix)
library(magrittr)
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)
ncpus <- 8

stimdat <- readRDS(file = "data/hvtn105_110subsets.rds")
stimdat$visitStimCell <- factor(stimdat$visitStimCell)
stimdat$stim <- factor(stimdat$stim)

# Analysis Settings ---------
configurations <- expand.grid(disp = c(10, 50),
                              iterations = c(30),
                              fseed = 1:20)
config <- configurations[setting, ]
iterations <- config[["iterations"]]
fseed <- config[["fseed"]]
disp <- config[["disp"]]
nPosteriors <- 4

control = flowReMix_control(updateLag = round(iterations / 3),
                            nPosteriors = nPosteriors,
                            maxDispersion = 10^3 * disp,
                            isingInit = -5,
                            intSampSize = 50,
                            ncores = ncpus,
                            seed = fseed,
                            isingWprior = FALSE,
                            isingStabilityReps = 100)

# stimdat <- subset(stimdat, visitStimCell %in% levels(visitStimCell)[221:240])
stimdat$ptid <- factor(stimdat$ptid)
stimdat$visitStimCell <- factor(stimdat$visitStimCell)
fit <- flowReMix(cbind(count, parentcount - count) ~ stim,
                subject_id = ptid,
                cell_type = visitStimCell,
                cluster_variable = stim,
                data = stimdat,
                covariance = "sparse",
                ising_model = "sparse",
                regression_method = "robust",
                cluster_assignment = TRUE,
                iterations = iterations,
                parallel = TRUE,
                verbose = TRUE,
                control = control)

fit$data <- NULL
filename <- paste("results/HVTN105_C_110subsets",
                  "TEST_iterations", iterations,
                  "disp", disp,
                  "nPost", nPosteriors,
                  "seed", fseed,
                  ".rds", sep = "")
saveRDS(fit, file = filename)

# plot(fit, type = "scatter", target = rx) +
#   theme(legend.position = "bottom") + theme_bw()
# plot(fit$isingStability)
# plot(fit, type = "boxplot", target = rx, jitter = TRUE,
#      normWeights = TRUE) +
#   theme(legend.position = "bottom")
#
# fit$data$targetI <- fit$data$rx == levels(fit$data$rx)[1]
# roctab <- summary(fit, type = "ROC", target = targetI)
# roctab[order(roctab$auc, decreasing = TRUE), ]
#
#
