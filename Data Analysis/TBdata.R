geomean <- function(x) exp(mean(log(x)))
library(flowReMix)
# Loading data -------------------------
tbdat <- readRDS("data/tb_rozot_booleans.RDS")
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
# tbdat$subset <- paste(tbdat$parent, tbdat$population, sep = "/")
# tbdat <- subset(tbdat, stim != "EBV")
# tbdat$stimtemp <- tbdat$stim
# tbdat$stim[tbdat$stim %in% c("P1", "P2", "P3")] <- "P"
# datlist <- list()
# stims <- unique(tbdat$stim)
# stims <- stims[stims != "UNS"]
# for(i in 1:length(stims)) {
#   temp <- subset(tbdat, stim %in% c(stims[i], "UNS"))
#   temp$stim[temp$stim == "UNS"] <- "ctrl"
#   temp$stimgroup <- stims[i]
#   datlist[[i]] <- temp
# }
# tbdat <- do.call("rbind", datlist)
# tbdat$subset <- paste(tbdat$stimgroup, tbdat$subset, sep = "/")
# nonzerocounts <- by(tbdat, tbdat$subset, function(x) mean(x$count > 0))
# nonzerocounts <- data.frame(names(nonzerocounts), as.numeric(unlist(nonzerocounts)))

# Defining subsets wo stim --------------
tbdat$subset <- paste(tbdat$parent, tbdat$population, sep = "/")
tbdat <- subset(tbdat, stim != "EBV")
tbdat$stim[tbdat$stim == "UNS"] <- "ctrl"
tbdat$stim <- factor(tbdat$stim, levels = c("ctrl", "MP", "Mtbaux", "P1", "P2", "P3"))
nonzerocounts <- by(tbdat, tbdat$subset, function(x) mean(x$count > 0))
nonzerocounts <- data.frame(names(nonzerocounts), as.numeric(unlist(nonzerocounts)))

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
library(ggplot2)
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


# Analysis ----------------------------------
npost <- 20
control <- flowReMix_control(updateLag = 15, nsamp = 50, initMHcoef = 1,
                             nPosteriors = npost, centerCovariance = TRUE,
                             maxDispersion = 10^3, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             lastSample = round(40 / npost), isingInit = -log(99),
                             initMethod = "sparse")

# tempdat$stim <- tempdat$stimtemp
# tempdat$stim[tempdat$stim == "UNS"] <- "aUNS"
# tempdat$stim <- factor(tempdat$stim, levels = sort(unique(tempdat$stim)))
# fit <- flowReMix(cbind(count, parentcount - count) ~ stim,
#                  subject_id = ptid,
#                  cell_type = subset,
#                  cluster_variable = stim,
#                  data = tempdat,
#                  covariance = "sparse",
#                  ising_model = "sparse",
#                  regression_method = "sparse",
#                  iterations = 25,
#                  parallel = TRUE,
#                  verbose = TRUE, control = control)
load(file = "data analysis/results/TBsep5.Robj")
# load(file = "data analysis/results/TBsep7robust.Robj")
# load(file = "data analysis/results/TBsep8robust.Robj")

# Scatter plots with posteriors ---------------
library(cowplot)
post <- fit$posteriors
outcome <- by(tempdat, tempdat$ptid, function(x) unique(x$type))
outcome <- data.frame(ptid = names(outcome), outcome = as.character(outcome))
post <- merge(post, outcome)
outcome <- post[, ncol(post)]
scatter <- plot(fit, type = "scatter", target = outcome,
                ncol = 10)
# save_plot(scatter, filename = "figures/TBscatter.pdf",
#           base_width = 20, base_height = 20)

# Fitting ROC curves -----------------
rocTable <- rocTable(fit, outcome, pvalue = "wilcoxon")
post <- fit$posteriors[, -1]
level <- 0.1
nresponders <- apply(post, 2, function(x) cummean(sort(1 - x)))
select <- nresponders[2, ] < level
rocTable$qvalue <- NA
rocTable$qvalue[select] <- p.adjust(rocTable$pvalue[select], method = "BH")
rocTable[order(rocTable$auc, decreasing = TRUE), ]
sum(rocTable$qvalue < 0.05, na.rm = TRUE)
sum(rocTable$qvalue < 0.1, na.rm = TRUE)

# Graph -----------------
load("data analysis/results/TBising5.Robj")
plot(stability, fill = rocTable$auc, threshold = 0.5)

# Boxplot by graph clusters -------------
groups <- getGraphComponents(stability, threshold = 0.5, minsize = 3)
plot(fit, type = "boxplot", target = outcome, groups = groups,
     test = "t-test", ncol = 2)

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
weights$Functionality <- rep(1, length(subsets))

group <- outcome
stimgroups <- sapply(subsets, function(x) strsplit(x, "/")[[1]][[1]])
stimnames <- unique(stimgroups)
stimgroups <- lapply(stimnames, function(x) subsets[stimgroups == x])
names(stimgroups) <- stimnames
plot(fit, type = "boxplot", weights = weights,
     target = group,
     test = "wilcoxon",
     groups = stimgroups)

cellgroups <- sapply(subsets, function(x) strsplit(x, "/")[[1]][[2]])
cellnames <- unique(cellgroups)
cellgroups <- lapply(cellnames, function(x) subsets[cellgroups == x])
names(cellgroups) <- cellnames
plot(fit, type = "boxplot", target = group, test = "wilcoxon",
     groups = cellgroups, ncol = 3, weights = weights)

stimcell <- sapply(subsets, function(x) paste(strsplit(x, "/")[[1]][1:2], collapse = "/"))
scnames <- unique(stimcell)
stimcell <- lapply(scnames, function(x) subsets[stimcell == x])
names(stimcell) <- scnames
stimcell <- stimcell[sapply(stimcell, function(x) length(x) > 0)]
plot(fit, type = "boxplot", target = group, test = "wilcoxon",
     groups = stimcell, ncol = 4)
