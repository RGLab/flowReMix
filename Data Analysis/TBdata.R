library(dplyr)
library(flowReMix)
library(ggplot2)

assign <- function(x) {
  x$prop <- x$count / x$parentcount
  assign <- as.numeric(by(x, x$subset, function(y) max(y$prop[y$stim != "aUNS"]) > min(y$prop[y$stim == "aUNS"])))
  assign[assign == 1] <- -1
  result <- data.frame(ptid = x$ptid[1], subset = unique(x$subset), assign = assign)
  return(result)
}

geomean <- function(x) exp(mean(log(x)))
# Loading data -------------------------
tbdat <- readRDS("data/TB_Rozot_booleans.rds")
names(tbdat) <- tolower(names(tbdat))
tbdat$ptid <- sapply(strsplit(tbdat$samplename, "_"), function(x) x[[1]])
tbdat$count[is.na(tbdat$count)] <- 0

# Parsing Expressions ----------------
expressions <- unique(tbdat$pop)
newexp <- character(length(expressions))
for(i in 1:length(expressions)) {
  e <- expressions[i]
  e <- strsplit(e, "&")[[1]]
  res <- c()
  for(j in 1:length(e)) {
    sube <- strsplit(e[j], "-")[[1]][1]
    if(substr(sube, 1, 1) != "!"){
      res <- c(res, sube)
    }
  }
  newexp[i] <- paste(paste(res, collapse = "+"), "+", sep = "")
}
map <- cbind(expressions, newexp)
for(i in 1:nrow(map)) {
  tbdat$population[tbdat$population == map[i, 1]] <- map[i, 2]
}

# Defining subsets w stim ---------------------
tbdat$subset <- paste(tbdat$parent, tbdat$population, sep = "/")
tbdat <- subset(tbdat, stim != "EBV")
tbdat$stimtemp <- tbdat$stim
tbdat$stim[tbdat$stim %in% c("P1", "P2", "P3")] <- "MP"
datlist <- list()
stims <- unique(tbdat$stim)
stims <- stims[stims != "UNS"]
for(i in 1:length(stims)) {
  temp <- subset(tbdat, stim %in% c(stims[i], "UNS"))
  temp$stim[temp$stim == "UNS"] <- "ctrl"
  temp$stimgroup <- stims[i]
  datlist[[i]] <- temp
}
tbdat <- do.call("rbind", datlist)
tbdat$subset <- paste(tbdat$stimgroup, tbdat$subset, sep = "/")
nonzerocounts <- by(tbdat, tbdat$subset, function(x) mean(x$count > 0))
nonzerocounts <- data.frame(names(nonzerocounts), as.numeric(unlist(nonzerocounts)))

# Defining subsets wo stim --------------
# tbdat$subset <- paste(tbdat$parent, tbdat$population, sep = "/")
# tbdat <- subset(tbdat, stim != "EBV")
# tbdat$stim[tbdat$stim == "UNS"] <- "ctrl"
# tbdat$stim <- factor(tbdat$stim, levels = c("ctrl", "MP", "Mtbaux", "P1", "P2", "P3"))
# nonzerocounts <- by(tbdat, tbdat$subset, function(x) mean(x$count > 0))
# nonzerocounts <- data.frame(names(nonzerocounts), as.numeric(unlist(nonzerocounts)))

# Choosing subset of data for analysis -----------------
countkeep <- nonzerocounts[nonzerocounts[, 2] >= 0.2, 1]
popkeep <- c("cd4", "MAIT", "NKrainbow", "DCs", "Bcells", "conCD*")
stimkeep <- c("P", "MP", "Mtbaux")
# tempdat <- subset(tbdat, (subset %in% countkeep) & (parent %in% popkeep) & (stimgroup %in% stimkeep))
tempdat <- subset(tbdat, (subset %in% countkeep) & (parent %in% popkeep))
#tempdat$stim <- factor(tempdat$stim, levels = c("ctrl", stimkeep))
tempdat$subset <- factor(tempdat$subset)
rm(tbdat)

# Visuaizations ------------------------------
tempdat$prop <- tempdat$count / tempdat$parentcount #+ 1 / tempdat$parentcount
subctrl <- subset(tempdat, stim == "ctrl")
substim <- subset(tempdat, stim != "ctrl")
wide <- merge(subctrl, substim, by = c("ptid", "subset"))
# pdf("figures/TBscat.pdf", width = 15, height = 15)
# ggplot(wide) +
#   geom_point(aes(x = prop.x, y = prop.y, col = type.x, shape = type.x)) +
#   facet_wrap(~ subset, scales = "free", nrow = floor(sqrt(length(unique(wide$subset))))) +
#   geom_abline(intercept = 0, slope = 1) +
#   theme_bw()
# dev.off()


set.seed(100)
seed = round(runif(50, min=1, max=500))
configurations = expand.grid(preAssignCoefs = 1,
                             isingWprior=FALSE,
                             markovChainEM=TRUE,
                             zeroPosteriorProbs = FALSE,
                             prior = 0,
                             sampleNew=FALSE,
                             seed=seed
                             )
args = as.numeric(commandArgs(trailingOnly = TRUE)[[2]])
settings = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
cpus = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
# args=1
# settings=1
# cpus=6
index = (args-1)*1 + settings
config = configurations[index,]
niter = 80
npost = 40


# Analysis ----------------------------------
control <- flowReMix_control(updateLag = 15, nsamp = 100, initMHcoef = 1,
                             nPosteriors = npost, centerCovariance = TRUE,
                             maxDispersion = 10^3, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             lastSample = NULL, isingInit = -log(99),
                             initMethod = "sparse",
                             preAssignCoefs = 1, prior = prior,
                             ncores = cpus,
                             keepEach = 1,
                             markovChainEM = TRUE,seed = seed,
                             isingWprior = FALSE,
                             zeroPosteriorProbs = FALSE,sampleNew = FALSE)

tempdat$stim <- tempdat$stimtemp
tempdat$stim[tempdat$stim == "UNS"] <- "aUNS"
tempdat$stim <- factor(tempdat$stim, levels = sort(unique(tempdat$stim)))
library(dplyr)
library(data.table)
# temptemp = tempdat%>%filter(stimgroup%like%"MP",subset%like%"NKrainbow")
# temptemp$subset = factor(temptemp$subset)
# temptemp$stim = factor(temptemp$stim)

fit <- flowReMix(cbind(count, parentcount - count) ~ stim,
                 subject_id = ptid,
                 cell_type = subset,
                 cluster_variable = stim,
                 data = as.data.frame(tempdat),
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "robust",
                 cluster_assignment = TRUE,
                 iterations = niter,
                 parallel = TRUE,
                 verbose = FALSE, control = control)

# add_ptid <- function(x, subject_id) {
#   x$subject_id <- match.call()$subject_id
#   return(x)
# }

# filenames <- as.list(dir(path = 'data analysis/results', pattern="TBdat3_perm_*"))
# keep <- which(apply(matrix(as.character(groups), nrow = 2), 2, function(x) any(sapply(x, function(y) grepl("MP.NKrainbow", y)))))
# filenum <- as.numeric(sapply(sapply(sapply(filenames, function(x) strsplit(x, "perm")), function(x) x[2]), function(x) strsplit(x, "rds")))
# filenames <- lapply(filenames, function(x) paste0('data analysis/results/', x))[filenum %in% keep][9]
filenames <- as.list(dir(path = 'data analysis/results', pattern="TBdat4_full_*"))
iter30 <- sapply(filenames, function(x) grepl("niter30", x))
post10 <-sapply(filenames, function(x) grepl("npost10", x))
filenames <- lapply(filenames, function(x) paste0('data analysis/results/', x))#[!iter30 & !post10]
post <- list()
for(i in 1:length(filenames)) {
  fit <- readRDS(file = filenames[[i]])
  post[[i]] <- fit$posteriors[, -1]
}
post <- Reduce("+", post) / length(filenames)
fit$posteriors[, -1] <- post

# stimcell  = lapply(split(tempdat$subset,interaction(factor(tempdat$stimgroup):factor(tempdat$parent))),unique)
# scboxplot <- plot(fit, type = "boxplot",
#                   target = type, test = "wilcoxon",
#                   #weights = weights,
#                   groups = stimcell, ncol = 4, jitter=TRUE)
# scboxplot
# scat <- plot(fit, type = "scatter", target = type)
# ggsave("figures/temptbscatter.pdf", plot = scat, width = 15, height = 15)


# Scatter plots with posteriors ---------------
library(cowplot)
post <- fit$posteriors
outcome <- by(tempdat, tempdat$ptid, function(x) unique(x$type))
outcome <- data.frame(ptid = names(outcome), outcome = as.character(outcome))
post <- merge(post, outcome)
outcome <- post[, ncol(post)]
# scatter <- plot(fit, type = "scatter", target = outcome,
#                 ncol = 10)
# save_plot(scatter, filename = "figures/TBscatter.pdf",
#           base_width = 20, base_height = 20)

# Fitting ROC curves -----------------
exclude <- sapply(fit$coefficients, function(x) max(x[-1]) < 0)
rocTable <- summary(fit, type = "ROC", test = "wilcoxon",
                    target = type)
post <- fit$posteriors[, -1]
level <- .50
nresponders <- apply(post, 2, function(x) cummean(sort(1 - x)))
select <- exclude == 0  & nresponders[1, ] < level
rocTable$qvalue <- NA
rocTable$qvalue[select] <- p.adjust(rocTable$pvalue[select], method = "BH")
rocTable[order(rocTable$auc, decreasing = TRUE), ]
sum(rocTable$qvalue < 0.05, na.rm = TRUE)
sum(rocTable$qvalue < 0.1, na.rm = TRUE)
# flowReMix:::summary.flowReMix(fit,target = outcome,type=c("ROC")) %>%
#   filter(subset %in% names(which(select))) %>%
#   arrange(-auc) %>% mutate(qvalue = p.adjust(pvalue,"BH")) %>% filter(responseProb>0.5,qvalue<0.1)

# Graph -----------------
stab <- fit$isingStability
plot(stab, threshold = 0.25, fill = rocTable$auc)

# Boxplot by graph clusters -------------
# groups <- getGraphComponents(stability, threshold = 0.5, minsize = 3)
# plot(fit, type = "boxplot", target = outcome, groups = groups,
#      test = "logistic", ncol = 2)

# Boxplots by categories -----------
subsets <- names(fit$posteriors)[-1]
functions <- sapply(subsets, function(x) strsplit(x, "/")[[1]][3])
functions <- unlist(lapply(functions, function(x) strsplit(x, "+", fixed = TRUE)[[1]]))
M <- length(unique(functions))
nfunctions <- sapply(subsets, function(x) strsplit(x, "/")[[1]][3])
nfunctions <- sapply(nfunctions, function(x) length(strsplit(x, "+", fixed = TRUE)[[1]]))
poly <- nfunctions / choose(rep(M, length(nfunctions)), nfunctions)
weights <- list()
# weights$polyfunctionality <- poly
# weights$Functionality <- rep(1, length(subsets))
# weights$weightedAvg <- apply(fit$posteriors[, -1], 2, sd)
weights <- list()
weights$exc <- rep(1, length(fit$coefficients))
#weights$exc[exclude] <- 0

allbox <- plot(fit, type = "boxplot",
                target = type, weights = weights,
                test = "wilcoxon",
                one_sided = TRUE,
                groups = "all", jitter = TRUE)
allbox
# save_plot(allbox, filename = "figures/TBallboxplot2.pdf",
#           base_height = 4, base_width = 8)

#
# group <- outcome
# stimgroups <- sapply(subsets, function(x) strsplit(x, "/")[[1]][[1]])
# stimnames <- unique(stimgroups)
# stimgroups <- lapply(stimnames, function(x) subsets[stimgroups == x])
# names(stimgroups) <- stimnames
stimgroups  = lapply(split(fit$data$subset,tempdat$stimgroup),unique)
stimbox <- plot(fit, type = "boxplot", weights = weights,
                target = type, test = "wilcoxon",
                # one_sided = TRUE,
                jitter = TRUE,
                groups = stimgroups)
stimbox
# save_plot(stimbox, filename = "figures/TBstimBoxplots2.pdf",
#           base_height = 4, base_width = 8)

# cellgroups <- sapply(subsets, function(x) strsplit(x, "/")[[1]][[2]])
# cellnames <- unique(cellgroups)
# cellgroups <- lapply(cellnames, function(x) subsets[cellgroups == x])
# names(cellgroups) <- cellnames
cellgroups  = lapply(split(tempdat$subset,tempdat$parent),unique)
cellbox <- plot(fit, type = "boxplot",
                target = type, test = "wilcoxon",
     groups = cellgroups, ncol = 3, weights = weights,
     jitter=TRUE)
cellbox
# save_plot(cellbox, filename = "figures/TBparentBoxplots2.pdf",
#           base_height = 4, base_width = 8)

stimcell <- sapply(subsets, function(x) paste(strsplit(x, "/")[[1]][1:2], collapse = "/"))
scnames <- unique(stimcell)
stimcell <- lapply(scnames, function(x) subsets[stimcell == x])
names(stimcell) <- scnames
stimcell <- stimcell[sapply(stimcell, function(x) length(x) > 0)]
stimcell  = lapply(split(tempdat$subset,interaction(factor(tempdat$parent):factor(tempdat$stimgroup))),unique)

scboxplot <- plot(fit, type = "boxplot",
                  target = type, test = "wilcoxon",
                  weights = weights,
                  groups = stimcell, ncol = 3, jitter=TRUE)
scboxplot
# save_plot(scboxplot, filename = "figures/TBstimParentBoxplot2.pdf",
#           base_height = 5, base_width = 10)
saveRDS(fit, file = paste0("TBfit_",index%%1+1,"_seed_",seed,".rds"))
system("sync")
