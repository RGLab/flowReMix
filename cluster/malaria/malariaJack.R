args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))

ncores <- 4
library(pryr)
library(flowReMix)
# Malaria dataset ----------------------------
load("data/malaria.rda")
populations <- unique(malaria$population)
parents <- unique(malaria$parent)
leaves <- populations[!(populations %in% parents) ]
malaria <- subset(malaria, population %in% leaves)
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

# Analysis -----------------------
library(flowReMix)
control <- flowReMix_control(updateLag = 15, keepEach = 5, nsamp = 50, initMHcoef = 2,
                             nPosteriors = 20, centerCovariance = TRUE,
                             maxDispersion = 500, minDispersion = 10^8,
                             randomAssignProb = 10^-8, intSampSize = 100,
                             lastSample = 10, isingInit = -log(95), ncores = ncores,
                             initMethod = "robust")

tempdat <- subset(malaria, TRUE)
tempdat$time <- tempdat$visitno
tempdat$subset <- factor(as.character(tempdat$subs))
tempdat$stim[tempdat$stim]
tempdat$trt <- 1
tempdat$visitno <- factor(as.character(tempdat$visitno), levels = c("Day 0", "Day 9", "pos",
                                                                    "Day 28", "Day 56", "Day 168"))
print(nrow(tempdat))
print(subject)
subject <- as.numeric(subject)
tempdat <- subset(tempdat, tempdat$ptid != unique(tempdat$ptid)[subject])
print(nrow(tempdat))

system.time(fit <- flowReMix(cbind(count, parentcount - count) ~
                               stim * visitno,
                             subject_id = ptid,
                             cell_type = subset,
                             cluster_variable = visitno,
                             data = tempdat,
                             covariance = "sparse",
                             ising_model = "sparse",
                             regression_method = "robust",
                             iterations = 25,
                             parallel = TRUE,
                             verbose = TRUE, control = control))
save(fit, file = paste("malariaJack2 subj", subject, ".Robj", sep = ""))

