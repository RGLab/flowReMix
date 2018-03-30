cpus <- 8
print(cpus)

args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)

getExpression <- function(str) {
  first <- substr(str, 1, 7)
  second <- substr(str, 8, nchar(str))
  second <- strsplit(second, "")[[1]]
  seperators <- c(0, which(second %in% c("-", "+")))
  expressed <- list()
  for(i in 2:length(seperators)) {
    if(second[seperators[i]] == "+") {
      expressed[[i]] <- paste(second[(seperators[(i - 1)] + 1) : seperators[i]], collapse = '')
    }
  }

  expressed <- paste(unlist(expressed), collapse = '')
  expressed <- paste(first, expressed, sep = '')
  return(expressed)
}

# Loading Data --------------------------------
# hvtn <- read.csv(file = "data/merged_505_stats.csv")
# names(hvtn) <- tolower(names(hvtn))
# hvtn <- subset(hvtn, !is.na(ptid))
# saveRDS(hvtn, file = "data/505_stats.rds")

# Getting marginals -----------------------------
library(flowReMix)
hvtn <- readRDS(file = "data/505_stats.rds")
length(unique(hvtn$name))
length(unique(hvtn$ptid))
length(unique(hvtn$population))
unique(hvtn$population)
unique(hvtn$stim)
nchars <- nchar(as.character(unique(hvtn$population)))
#marginals <- unique(hvtn$population)[nchars < 26]
marginals <- unique(hvtn$population)[nchars == 26]
marginals <- subset(hvtn, population %in% marginals)
marginals <- subset(marginals, stim %in% c("negctrl", "VRC ENV A",
                                           "VRC ENV B", "VRC ENV C",
                                           "VRC GAG B", "VRC NEF B",
                                           "VRC POL 1 B", "VRC POL 2 B"))
marginals <- subset(marginals, !(population %in% c("4+", "8+")))
marginals <- subset(marginals, !(population %in% c("8+/107a-154-IFNg-IL2-TNFa-", "4+/107a-154-IFNg-IL2-TNFa-")))
marginals$stim <- factor(as.character(marginals$stim))
marginals$population <- factor(as.character(marginals$population))

# Descriptives -------------------------------------
library(ggplot2)
marginals$prop <- marginals$count / marginals$parentcount
# ggplot(marginals) + geom_boxplot(aes(x = population, y = log(prop), col = stim))

require(dplyr)
negctrl <- subset(marginals, stim == "negctrl")
negctrl <- summarize(group_by(negctrl, ptid, population), negprop = mean(prop))
negctrl <- as.data.frame(negctrl)
marginals <- merge(marginals, negctrl, all.x = TRUE)

# ggplot(subset(marginals, stim != "negctrl" & parent == "4+")) +
#   geom_point(aes(x = log(negprop), y = log(prop)), size = 0.25) +
#   facet_grid(stim ~ population, scales = "free") +
#   theme_bw() +
#   geom_abline(intercept = 0, slope = 1)

# Setting up data for analysis ---------------------------
unique(marginals$stim)
gag <- subset(marginals, stim %in% c("VRC GAG B", "negctrl"))
gag$subset <- factor(paste("gag", gag$population, sep = "/"))
gag$stimGroup <- "gag"
pol <-subset(marginals, stim %in% c("negctrl", "VRC POL 1 B", "VRC POL 2 B"))
pol$subset <- factor(paste("pol", pol$population, sep = "/"))
pol$stimGroup <- "pol"
env <- subset(marginals, stim %in% c("negctrl", "VRC ENV C", "VRC ENV B", "VRC ENV A"))
env$subset <- factor(paste("env", env$population, sep = "/"))
env$stimGroup <- "env"
nef <- subset(marginals, stim %in% c("negctrl", "VRC NEF B"))
nef$subset <- factor(paste("nef", nef$population, sep = "/"))
nef$stimGroup <- "nef"
subsetDat <- rbind(gag, pol, env, nef)
subsetDat$stim <- as.character(subsetDat$stim)
subsetDat$stim[subsetDat$stim == "negctrl"] <- 0
subsetDat$stim <- factor(subsetDat$stim)

# Converting subset names ------------------
subsets <- as.character(unique(subsetDat$subset))
expressed <- sapply(subsets, getExpression)
map <- cbind(subsets, expressed)
subsetDat$subset <- as.character(subsetDat$subset)
for(i in 1:nrow(map)) {
  subsetDat$subset[which(subsetDat$subset == map[i, 1])] <- map[i, 2]
}
subsetDat$subset <- factor(subsetDat$subset)

# Getting outcomes -------------------------------
# treatmentdat <- read.csv(file = "data/rx_v2.csv")
# names(treatmentdat) <- tolower(names(treatmentdat))
# treatmentdat$ptid <- factor(gsub("-", "", (treatmentdat$ptid)))
# treatmentdat <- subset(treatmentdat, ptid %in% unique(subsetDat$ptid))

# Finding problematic subsets?
keep <- by(subsetDat, list(subsetDat$subset), function(x) mean(x$count > 1) > 0.02)
# keep <- by(subsetDat, list(subsetDat$subset), function(x) sum(x$count > 5) > 2)
keep <- names(keep[sapply(keep, function(x) x)])
#result$subsets[result$qvals < 0.1] %in% keep
subsetDat <- subset(subsetDat, subset %in% keep)
subsetDat$subset <- factor(as.character(subsetDat$subset))

configurations <- expand.grid(method = c("MC"),
                              seed = 1:50,
                              prior = c(0),
                              niter = c(60),
                              includeBatch = FALSE)
config <- configurations[setting, ]
print(config)
niter <- config[["niter"]]
seed <- config[["seed"]]
prior <- config[["prior"]]
method <- config[["method"]]
includeBatch <- config[["includeBatch"]]
if(method == "MC") {
  npost <- 1 #3
  lag <- round(niter / 3)
  keepeach <- 5
  mcEM <- TRUE
}

# Keeping max collection -----
subsetDat$batch <- factor(subsetDat$batch..)
subsetDat$stimGroup <- factor(subsetDat$stimGroup)
subsetDat <- data.frame(subsetDat %>% group_by(ptid,population,stim,stimGroup,parent) %>% filter(collection.num==max(collection.num)))
subsetDat$batch <- factor(as.character(subsetDat$batch), levels = unique(as.character(subsetDat$batch)))

# Merging with hla -------
hla <- read.csv("data/p505_HLA_category.csv")
names(hla) <- tolower(names(hla))
subsetDat <- merge(subsetDat, hla, all.x = TRUE, all.y = FALSE)
# hla_status <- subsetDat %>% select(ptid, hla_cat) %>% unique() %>%
#   select(hla_cat) %>% table(useNA = "ifany")
subsetDat$hla_cat <- as.character(subsetDat$hla_cat)
subsetDat$hla_cat[is.na(subsetDat$hla_cat)] <- "miss"
subsetDat$hla_cat <- as.factor(subsetDat$hla_cat)

# Fitting the model ------------------------------
library(flowReMix)
control <- flowReMix_control(updateLag = lag, nsamp = 50,
                             keepEach = keepeach, initMHcoef = 2.5,
                             nPosteriors = npost, centerCovariance = FALSE,
                             maxDispersion = 10^4 * 5, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             seed = seed, zeroPosteriorProbs = FALSE,
                             ncores = cpus, preAssignCoefs = 1,
                             prior = prior, isingWprior = FALSE,
                             markovChainEM = mcEM,
                             initMethod = "robust",
                             learningRate = 0.6, keepWeightPercent = 0.9)

fit <- flowReMix(cbind(count, parentcount - count) ~ stim * hla_cat,
                 subject_id = ptid,
                 cell_type = subset,
                 cluster_variable = stim,
                 data = subsetDat,
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "robust",
                 iterations = niter,
                 parallel = TRUE, keepSamples = TRUE,
                 cluster_assignment = TRUE,
                 verbose = TRUE, control = control)

file <- paste("results/hvtn_hla_1_1045_niter", niter, "npost", npost, "seed", seed, "prior", prior, method, ".rds", sep = "")
print(file)
saveRDS(object = fit, file = file)
stab <- stabilityGraph(fit, type = "ising", cpus = cpus, AND = TRUE,
                       gamma = 0.25, reps = 50, cv = FALSE)
fit$stabilityGraph <- stab
fit$randomEffectSamp <- NULL
fit$assignmentList <- NULL
fit$data <- NULL
saveRDS(object = fit, file = file)

