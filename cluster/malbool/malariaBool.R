library(flowReMix)
library(magrittr)
library(dplyr)
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)
ncores <- 6

# Loading data ----
malbool <- readRDS("data/malaria_booleans.rds")
names(malbool) <- tolower(names(malbool))

# Setting up stimulations -----
malbool$subset <- with(malbool, (interaction(parent, population, sep = "/")))
malbool$stim <- as.character(malbool$stim)
malbool$stim[malbool$stim == "Cells Only"] <- "ctrl"
malbool$stim[malbool$stim == "ctr SPZ"] <- "ctrl"
malbool$stim[malbool$stim == "PfSPZ"] <- "SPZ"
malbool$stimgroup <- "SPZ"
malbool$stimgroup[malbool$stim %in% c("uRBC", "PfRBC")] <- "RBC"
malbool$stim[malbool$stim == "uRBC"] <- "ctrl"
malbool$stim <- factor(malbool$stim)

# Defining and screening subsets ---
malbool$subset <- factor(with(malbool, (interaction(stimgroup, subset, sep = "/"))))
counts <- by(malbool, list(malbool$subset), function(x) x$count)
dropsubsets <- sapply(counts, function(x) mean(x >= 2) < 0.15)
names(counts)[!dropsubsets]
malbool$subset <- as.character(malbool$subset)
malbool <- subset(malbool, subset %in% names(counts)[!dropsubsets])
malbool$subset <- factor(malbool$subset)

# Setting correct order for visits ---
malbool$visitno <- as.character(malbool$visitno)
malbool$visitno[malbool$visitno %in% c("Day 13-Pos (day of blood stage parasitemia)",
                                       "Day 14-Pos (day of blood stage parasitemia)",
                                       "Day 19-Pos (day of blood stage parasitemia)",
                                       "Day 11-Pos (day of blood stage parasitemia)",
                                       "Day 17-Pos (day of blood stage parasitemia)")] <- "pos"
malbool$visitno <- factor(malbool$visitno, levels = c("Day 0", "Day 9", "pos", "Day 28", "Day 56", "Day 168"))

# Analysis Parameters -----
configurations <- expand.grid(iterations = c(60),
                              mcEM = TRUE,
                              disp = c(10, 50),
                              npost = c(20),
                              seed = 1:30)
config <- configurations[setting, ]
niter <- config[["iterations"]]
lag <- round(niter / 3)
mcEM <- config[["mcEM"]]
npost <- config[["npost"]]
seed <- config[["seed"]]
disp <- config[["disp"]]
stabilityReps <- 50

# Control Object --------
control <- flowReMix_control(updateLag = lag, nsamp = 50,
                             keepEach = keepeach, initMHcoef = 2.5,
                             nPosteriors = npost, centerCovariance = FALSE,
                             maxDispersion = maxdisp * 1000, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             seed = seed, zeroPosteriorProbs = FALSE,
                             ncores = cpus, preAssignCoefs = 1,
                             prior = prior, isingWprior = FALSE,
                             markovChainEM = mcEM,
                             initMethod = "robust")

malbool$allstims <- malbool$subset %>% as.character() %>%
  strsplit("/") %>% sapply(function(x) paste(x[-1], collapse = "/")) %>%
  factor()

# Analysis --------
fit <- flowReMix(cbind(count, parentcount - count) ~ visitno * stim,
                 subject_id = ptid,
                 cell_type = allstims,
                 cluster_variable = visitno,
                 data = malbool,
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "robust",
                 cluster_assignment = TRUE,
                 iterations = niter,
                 parallel = TRUE,
                 verbose = TRUE, control = control)

file <- paste("results/malbool_joint",
              "seed", seed,
              "npost", npost,
              "niter", niter,
              ".rds", sep ="")
saveRDS(fit, file = file)











