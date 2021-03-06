library(flowReMix)
library(magrittr)
require(dplyr)
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)
cpus <- 2
print(cpus)

# Helper functions
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

getExpressionB <- function(str) {
  split <- strsplit(str, "/")
  parent <- split[[1]][1]
  pop <- split[[1]][2]
  split <- strsplit(pop, "+", fixed = TRUE) %>% unlist()
  for(i in 1:length(split)) {
    s <- split[i]
    if(i == length(split)) {
      last <- substr(s, nchar(s), nchar(s))
      if(last == "-") {
        split[i] <- ""
      } else {
        s <- strsplit(s, "-", fixed = TRUE) %>% unlist()
        split[i] <- s[length(s)]
      }
    } else {
      s <- strsplit(s, "-", fixed = TRUE) %>% unlist()
      split[i] <- s[length(s)]
    }
  }
  pop <- paste(split, collapse = "+")
  result <- paste(parent, pop, sep = "/")
  return(result)
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
                                           # "Empty Ad5 (VRC)"))
marginals <- subset(marginals, !(population %in% c("4+", "8+")))
marginals <- subset(marginals, !(population %in% c("8+/107a-154-IFNg-IL2-TNFa-", "4+/107a-154-IFNg-IL2-TNFa-")))
marginals$stim <- factor(as.character(marginals$stim))
marginals$population <- factor(as.character(marginals$population))

# Descriptives -------------------------------------
marginals$prop <- marginals$count / marginals$parentcount
# ggplot(marginals) + geom_boxplot(aes(x = population, y = log(prop), col = stim))

negctrl <- subset(marginals, stim == "negctrl")
negctrl <- summarize(group_by(negctrl, ptid, population), negprop = mean(prop))
negctrl <- as.data.frame(negctrl)
marginals <- merge(marginals, negctrl, all.x = TRUE)

# Converting subset names ------------------
subsets <- as.character(unique(marginals$population))
expressed <- sapply(subsets, getExpressionB)
map <- cbind(subsets, expressed)
marginals$population <- as.character(marginals$population)
for(i in 1:nrow(map)) {
  marginals$population[which(marginals$population == map[i, 1])] <- map[i, 2]
}
marginals$population <- factor(marginals$population)
marginals <- marginals %>% group_by(ptid, population, stim, parent) %>%
  filter(collection.num==max(collection.num)) %>% data.frame()

# Setting up data for analysis ---------------------------
# # Screening subsets based on mixed effect model ---------
# library(lme4)
# screenResults <- data.frame(subset = levels(marginals$population),
#                             pval = 1)
# for(i in 1:length(levels(marginals$population))) {
#   subdat <- subset(marginals, marginals$population == levels(marginals$population)[i])
#   subdat$od <- 1:nrow(subdat)
#   lmeFit <- NULL
#   try(lmeFit <- glmer(cbind(count, parentcount - count) ~ stim + (1|od) + (1|ptid),
#                   data = subdat, family = "binomial") %>% summary())
#   if(is.null(lmeFit)) {
#     next
#   }
#   nullFit <- NULL
#   try(nullFit <- glmer(cbind(count, parentcount - count) ~ (1|od) + (1|ptid),
#                        data = subdat, family = "binomial") %>% summary())
#   if(is.null(nullFit)) {
#     next
#   }
#
#   likRatio <- 2 * (lmeFit$logLik - nullFit$logLik)
#   lrDF <- nrow(coef(lmeFit)) - nrow(coef(nullFit))
#   likRatioPval <- pchisq(likRatio, lrDF, lower.tail = FALSE) %>% as.numeric()
#   screenResults[i, 2] <- likRatioPval
#   print(screenResults[i, ])
# }
# saveRDS(screenResults, file = "data/HVTN505jointScreen.rds")
screenResults <- readRDS(file = "data/HVTN505jointScreen.rds")
screenResults$qval <- p.adjust(screenResults$pval, method = "BH")
screenResults <- screenResults[order(screenResults$pval), ]

# Finding problematic subsets?
keep <- subset(screenResults, qval < 0.2)$subset
marginals <- subset(marginals, population %in% keep)
marginals$population <- factor(as.character(marginals$population))

# keep <- by(marginals, list(marginals$population), function(x) mean(x$count > 1) > 0.02)
# keep <- names(keep[sapply(keep, function(x) x)])
# marginals <- subset(marginals, population %in% keep)
# marginals$population <- factor(as.character(marginals$population))

configurations <- expand.grid(method = c("MC"),
                              seed = 1:50,
                              maxdisp = c(10, 50),
                              niter = c(60),
                              includeBatch = FALSE)
config <- configurations[setting, ]
print(config)
niter <- config[["niter"]]
seed <- config[["seed"]]
prior <- 0
maxdisp <- config[["maxdisp"]]
method <- config[["method"]]
includeBatch <- config[["includeBatch"]]
if(method == "MC") {
  npost <- 1
  lag <- round(niter / 3)
  keepeach <- 5
  mcEM <- TRUE
}

# Fitting the model ------------------------------
control <- flowReMix_control(updateLag = lag, nsamp = 50,
                             keepEach = keepeach, initMHcoef = 2.5,
                             nPosteriors = npost, centerCovariance = FALSE,
                             maxDispersion = maxdisp * 1000, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             seed = seed, zeroPosteriorProbs = FALSE,
                             ncores = cpus, preAssignCoefs = 1,
                             prior = prior, isingWprior = FALSE,
                             markovChainEM = mcEM,
                             initMethod = "robust",
                             learningRate = 0.6, keepWeightPercent = 0.9)

fit <- flowReMix(cbind(count, parentcount - count) ~ stim,
                 subject_id = ptid,
                 cell_type = population,
                 cluster_variable = stim,
                 data = marginals,
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "robust",
                 iterations = niter,
                 parallel = TRUE, keepSamples = TRUE,
                 cluster_assignment = TRUE,
                 verbose = TRUE, control = control)

file <- paste("results/hvtn_jointF_A",
              "_maxdisp", maxdisp,
              "_niter", niter,
              "npost", npost,
              "seed", seed,
              "prior", prior,
              method, ".rds", sep = "")
print(file)
saveRDS(object = fit, file = file)
stab <- stabilityGraph(fit, type = "ising", cpus = cpus, AND = TRUE,
                       gamma = 0.25, reps = 100, cv = FALSE)
fit$stabilityGraph <- stab
fit$randomEffectSamp <- NULL
fit$assignmentList <- NULL
fit$data <- NULL
saveRDS(object = fit, file = file)

