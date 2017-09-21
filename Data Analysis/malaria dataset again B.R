library(flowReMix)
# Malaria dataset ----------------------------
data(malaria)
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
malaria$stim <- factor(malaria$stim, levels = c("control", "PfSPZ", "PfRBC"))
isCytokine <- substring(malaria$population, nchar(malaria$population)) == "+"
malaria <- subset(malaria, isCytokine)
malaria$subset <- paste(malaria$visitno, "/", malaria$stimgroup, "/", malaria$population, sep = "")
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
control <- flowReMix_control(updateLag = 4, nsamp = 100, initMHcoef = 1,
                             nPosteriors = 2, centerCovariance = TRUE,
                             maxDispersion = 10^4, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             lastSample = 40, isingInit = -log(89),
                             initMethod = "robust")

tempdat <- subset(malaria, parent %in% c("4+"))
tempdat$subset <- factor(as.character(tempdat$subs))
tempdat$trt <- 1
system.time(fit <- flowReMix(cbind(count, parentcount - count) ~ trt * stim,
                 subject_id = ptid,
                 cell_type = subset,
                 cluster_variable = stim,
                 data = tempdat,
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "robust",
                 iterations = 8,
                 parallel = FALSE,
                 verbose = TRUE, control = control))
# save(fit, file = "data analysis/results/new malaria 4+ 4+CXCR5+.Robj")
load(file = "data analysis/results/malaria2_5_npost20niter20maxDisp1000.Robj")
load(file = "data analysis/results/malaria2_4_npost10niter20maxDisp100.Robj")

# ROC results ------------------------
infect <- by(malaria, malaria$ptid, function(x) unique(x$infection))
infect <- data.frame(ptid = names(infect), status = as.vector(infect))

plot(fit, type = "boxplot", target = infect[, 2], groups = "all",
     test = "wilcoxon")
rocs <- summary(fit, type = "ROC", target = )

