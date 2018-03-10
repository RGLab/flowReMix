library(flowReMix)
library(magrittr)
library(reshape2)
library(dplyr)

args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)
ncpus <- 4

filtered_data <- readRDS(file="data/hvtn105_filtered_data.rds")

# Defining (time X stim X cell) subsets
names(filtered_data) <- tolower(names(filtered_data))
unique(filtered_data$stim)
filtered_data <- filtered_data[, -14]
filtered_data$subset <- gsub(".", "", filtered_data$subset, fixed = "TRUE")
stimdat <- stimulationModel(filtered_data,
                            cell_type = subset,
                            stim_var = stim,
                            controls = "negctrl",
                            stim_groups = list(gag = "Gag-ZM96",
                                               vacenv = "92TH023-ENV",
                                               env = c("Env-1-ZM96", "Env-2-ZM96")))
unique(stimdat$stimCellType)
stimdat$visitStimCell <- interaction(stimdat$visitno, stimdat$stimCellType, sep = "/", drop = TRUE)

# Screening samples with very few cells ----
stimdat$totalCells <- stimdat$name %>% strsplit(split = "_") %>%
  sapply(function(x) x[2]) %>% as.numeric()
stimdat <- subset(stimdat, totalCells > 10^5)

# Screening samples with very small parent count ----
stimdat <- subset(stimdat, parentcount >= 100)

# Screening samples that are not paired ------
paired <- group_by(stimdat, visitStimCell, ptid) %>%
  summarize(pairedInd = sum(stim == "ctrl") > 0 & sum(stim != "ctrl") > 0)
stimdat <- merge(stimdat, paired, all.x = TRUE)
stimdat <- subset(stimdat, pairedInd)

countSummary <- group_by(stimdat, stimCellType, visitno, ptid) %>%
  summarize(nonzero = any(count > 0)) %>%
  group_by(stimCellType, visitno) %>% summarize(nonzero = sum(nonzero)) %>%
  group_by(stimCellType) %>% summarize(nonzero = mean(nonzero >= 10)) %>%
  data.frame() %>% subset(nonzero == 1)

stimdat <- subset(stimdat, stimCellType %in% countSummary$stimCellType)
stimdat$visitStimCell <- factor(as.character(stimdat$visitStimCell))
stimdat$stimCellType <- factor(stimdat$stimCellType)
levels(stimdat$visitStimCell)
stimdat <- data.frame(stimdat)

# Analysis Settings ---------
configurations <- expand.grid(pargroup = 1:13,
                              iterations = 40,
                              fseed = 1:4)
config <- configurations[setting, ]
iterations <- config[["iterations"]]
fseed <- config[["fseed"]]
pargroup <- config[["pargroup"]]

stimdat$par <- as.character(stimdat$visitStimCell) %>% strsplit("/") %>%
  sapply(function(x) x[[3]])
parentGroups <- unique(stimdat$par)
stimdat <- subset(stimdat, par == parentGroups[pargroup])
stimdat$visitStimCell <- factor(stimdat$visitStimCell)
stimdat$stim <- factor(stimdat$stim)

nPosteriors <- ceiling(length(levels(stimdat$visitStimCell)) / length(unique(stimdat$ptid))) + 1

control = flowReMix_control(updateLag = round(iterations / 3),
                            nPosteriors = nPosteriors,
                            maxDispersion = 10^4,
                            isingInit = -6,
                            intSampSize = 50,
                            ncores = ncpus,
                            seed = fseed,
                            isingWprior = FALSE,
                            isingStabilityReps = 100,
                            threads = ncpus * 2)

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

filename <- paste("results/HVTN105_A_parent_", parentGroups[pargroup],
                  "_iterations", iterations,
                  "nPost", nPosteriors,
                  "seed", fseed,
                  ".rds", sep = "")
saveRDS(fit, file = filename)
