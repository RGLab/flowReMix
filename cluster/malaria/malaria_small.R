library(flowReMix)
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)
cpus <- 1

# Malaria dataset ----------------------------
load("data/malaria.rda")
names(malaria)
table(malaria$experiment)
unique(malaria$ptid)
unique(malaria$population)
populations <- unique(malaria$population)
parents <- unique(malaria$parent)
leaves <- populations[!(populations %in% parents) ]
malaria <- subset(malaria, population %in% leaves)
unique(malaria$stim)
malaria$stimgroup[malaria$stim %in% c("PfRBC", "uRBC")] <- "RBC"
malaria$stimgroup[!(malaria$stim %in% c("PfRBC", "uRBC"))] <- "SPZ"
malaria$stim[malaria$stim == "uRBC"] <- "control"
malaria$stim[malaria$stim != "control"] <- "stim"
malaria$stim <- factor(malaria$stim, levels = c("control", "stim"))
isCytokine <- substring(malaria$population, nchar(malaria$population)) == "+"
malaria <- subset(malaria, isCytokine)
# malaria$subset <- paste(malaria$visitno, "/", malaria$stimgroup, "/", malaria$population, sep = "")
malaria$subset <- paste(malaria$stimgroup, "/", malaria$population, sep = "")
malaria$visitno <- factor(malaria$visitno)
malaria <- subset(malaria, parent %in% c("4+", "4+/CXCR5+", "4+/CXCR5+/PD-1+", "8+", "8+/CXCR5+"))

malaria$infection <- TRUE
malaria$infection[malaria$ptid %in% c("60061", "50071", "20003")] <- FALSE

# Screening low counts -------------------
countlist <- by(malaria, malaria$subset, function(x) x$count)
toRemove <- sapply(countlist, function(x) mean(x > 4) < 0.05)
toRemove <- names(countlist)[toRemove]
malaria <- subset(malaria, !(subset %in% toRemove))
malaria$subset <- factor(malaria$subset)
malaria$ptid <- factor(malaria$ptid)
malaria <- malaria[order(malaria$ptid, malaria$stimgroup), ]

# Further narrowing down to RBC ------
isRBC <- grepl("RBC", malaria$subset)
malaria <- subset(malaria, isRBC)
malaria$subset <- factor(malaria$subset)


# Analysis Settings -----
configurations <- expand.grid(method = c("MC"),
                              seed = 1:25,
                              maxdisp = c(10, 50),
                              niter = c(60))

config <- configurations[setting, ]
print(config)
niter <- config[["niter"]]
seed <- config[["seed"]]
prior <- 0
maxdisp <- config[["maxdisp"]]
method <- config[["method"]]
if(method == "MC") {
  npost <- 5
  lag <- round(niter / 3)
  keepeach <- 5
  mcEM <- TRUE
}


# Analysis -----------------------
library(flowReMix)
control <- flowReMix_control(updateLag = lag, nsamp = 50,
                             keepEach = keepeach, initMHcoef = 0.25,
                             nPosteriors = npost, centerCovariance = FALSE,
                             maxDispersion = maxdisp * 1000, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             seed = seed, zeroPosteriorProbs = FALSE,
                             ncores = cpus, preAssignCoefs = 1,
                             prior = prior, isingWprior = FALSE,
                             markovChainEM = mcEM,
                             initMethod = "robust")

tempdat <- subset(malaria, TRUE)
tempdat$time <- tempdat$visitno
tempdat$subset <- factor(as.character(tempdat$subs))
tempdat$stim[tempdat$stim]
tempdat$trt <- 1
tempdat$visitno <- as.character(tempdat$visitno)
tempdat <- subset(tempdat, visitno %in% c("Day 0", "Day 9", "Day 28", "Day 168"))
tempdat$visitno <- factor(as.character(tempdat$visitno), levels = c("Day 0", "Day 9", "Day 28", "Day 168"))
# saveRDS(tempdat, file = "data/small_malaria_marginal.rds")
system.time(fit <- flowReMix(cbind(count, parentcount - count) ~ stim * visitno,
                             subject_id = ptid,
                             cell_type = subset,
                             cluster_variable = visitno,
                             data = tempdat,
                             covariance = "sparse",
                             ising_model = "sparse",
                             regression_method = "robust",
                             iterations = niter,
                             parallel = TRUE,
                             verbose = TRUE, control = control))

file <- paste("results/malaria_marginal_small_RBC4+8+_A",
              "_maxdisp", maxdisp,
              "niter", niter,
              "npost", npost,
              "seed", seed,
              "prior", prior,
              method, ".rds", sep = "")
print(file)
fit$randomEffectSamp <- NULL
fit$assignmentList <- NULL
fit$data <- NULL
saveRDS(object = fit, file = file)
