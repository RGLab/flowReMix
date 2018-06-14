library(flowReMix)
library(magrittr)
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)
ncpus <- 4

stimdat <- readRDS(file = "data/hvtn105temp67subsets.rds")
stimdat$visitStimCell <- factor(stimdat$visitStimCell)
stimdat$stimCellType <- factor(stimdat$stimCellType)
stimdat$stim <- factor(stimdat$stim)

# Analysis Settings ---------
configurations <- expand.grid(disp = c(10, 50),
                              iterations = c(60),
                              fseed = 1:30)
config <- configurations[setting, ]
iterations <- config[["iterations"]]
fseed <- config[["fseed"]]
disp <- config[["disp"]]
nPosteriors <- 1

control = flowReMix_control(updateLag = round(iterations / 3),
                            nPosteriors = nPosteriors,
                            maxDispersion = 10^3 * disp,
                            isingInit = -4,
                            intSampSize = 50,
                            ncores = ncpus,
                            seed = fseed,
                            isingWprior = FALSE,
                            isingStabilityReps = 100)

# stimdat <- subset(stimdat, visitStimCell %in% levels(visitStimCell)[221:240])
stimdat$ptid <- factor(stimdat$ptid)
stimdat$stimVisit <- interaction(stimdat$stim, stimdat$visitno)
stimdat$stimVisit <- as.character(stimdat$stimVisit)
stimdat$stimVisit[stimdat$stimVisit %in% c("ctrl.5", "ctrl.7", "ctrl.9")] <- "ctrl"
stimdat$stimVisit <- factor(stimdat$stimVisit)
whichCtrl <- which(levels(stimdat$stimVisit) == "ctrl")
levels(stimdat$stimVisit) <- c("ctrl", levels(stimdat$stimVisit)[-whichCtrl]) %>% unique()
fit <- flowReMix(cbind(count, parentcount - count) ~ stimVisit,
                subject_id = ptid,
                cell_type = stimCellType,
                cluster_variable = stimVisit,
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
filename <- paste("results/HVTN105temp_stimVisitResponse_A",
                  "_iterations", iterations,
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
