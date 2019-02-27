library(flowReMix)
library(magrittr)
library(dplyr)
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)
ncores <- 1

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

# Screening based on mixed models ---------
# library(lme4)
# subsets <- levels(malbool$subset)
# screenResults <- data.frame(subset = subsets, pvalue = 1)
# for(i in 1:length(subsets)) {
#   pop <- subsets[[i]]
#   subdat <- subset(malbool, subset == pop)
#   subdat$sampind <- 1:nrow(subdat)
#
#   fullmod <- NULL
#   try(fullmod <- glmer(cbind(count, parentcount - count) ~ stim * visitno + (1|sampind) + (1|ptid),
#                    family = "binomial",
#                    data = subdat))
#   if(is.null(fullmod)) next
#   fullLogLik <- summary(fullmod)$logLik[[1]]
#
#   nullmod <- NULL
#   try(nullmod <- glmer(cbind(count, parentcount - count) ~ stim + (1|sampind) + (1|ptid),
#                    family = "binomial",
#                    data = subdat))
#   if(is.null(nullmod)) next
#   nullLoglik <- summary(nullmod)$logLik[[1]]
#
#   df <- nrow(summary(fullmod)$coefficients) - nrow(summary(nullmod)$coefficients)
#   chisqStat <- 2 * (fullLogLik - nullLoglik)
#   screenResults$pvalue[i] <- 1 - pchisq(chisqStat, df)
#   print(screenResults[i, ])
# }
# saveRDS(screenResults, file = "data/malboolScreen.rds")
screenResults <- readRDS(file = "data/malboolScreen.rds")
screenResults <- screenResults[order(screenResults$pvalue), ]
malbool <- subset(malbool, subset %in% screenResults[1:50, 1])
malbool$subset <- factor(malbool$subset)

# counts <- by(malbool, list(malbool$subset), function(x) x$count)
# dropsubsets <- sapply(counts, function(x) mean(x >= 2) < 0.15)
# names(counts)[!dropsubsets]
# malbool$subset <- as.character(malbool$subset)
# malbool <- subset(malbool, subset %in% names(counts)[!dropsubsets])

# Setting correct order for visits ---
malbool$visitno <- as.character(malbool$visitno)
malbool$visitno[malbool$visitno %in% c("Day 13-Pos (day of blood stage parasitemia)",
                                       "Day 14-Pos (day of blood stage parasitemia)",
                                       "Day 19-Pos (day of blood stage parasitemia)",
                                       "Day 11-Pos (day of blood stage parasitemia)",
                                       "Day 17-Pos (day of blood stage parasitemia)")] <- "pos"
malbool$visitno <- factor(malbool$visitno, levels = c("Day 0", "Day 9", "pos", "Day 28", "Day 56", "Day 168"))

# Analysis Parameters -----
configurations <- expand.grid(iterations = c(50),
                              mcEM = TRUE,
                              disp = c(50),
                              npost = c(8),
                              seed = 1:40,
                              visitno =  c("pos", "Day 56", "Day 168"))
config <- configurations[setting, ]
niter <- config[["iterations"]]
lag <- round(niter / 3)
mcEM <- config[["mcEM"]]
npost <- config[["npost"]]
seed <- config[["seed"]]
disp <- config[["disp"]]
visit = config[["visitno"]]
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
                             sampleNew = FALSE)

malbool$allstims <- malbool$subset %>% as.character() %>%
  strsplit("/") %>% sapply(function(x) paste(x[-1], collapse = "/")) %>%
  factor()

# Analysis --------
malbool <- subset(malbool, subset %in% levels(subset))
malbool$subset <- factor(malbool$subset)
malbool <- subset(malbool, as.character(visitno) == as.character(visit))
fit <- flowReMix(cbind(count, parentcount - count) ~ stim,
                 subject_id = ptid,
                 cell_type = subset,
                 cluster_variable = stim,
                 data = malbool,
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "robust",
                 cluster_assignment = TRUE,
                 iterations = niter,
                 parallel = TRUE,
                 verbose = TRUE, control = control)

file <- paste("results/malbool_screenB_",
              gsub(" ", "", visit),
              "_disp", disp,
              "seed", seed,
              "npost", npost,
              "niter", niter,
              ".rds", sep ="")
saveRDS(fit, file = file)











