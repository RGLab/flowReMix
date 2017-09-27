args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)
ncores <- 8

assign <- function(x) {
  x$prop <- x$count / x$parentcount
  assign <- as.numeric(by(x, x$subset, function(y) max(y$prop[y$stim != "aUNS"]) > min(y$prop[y$stim == "aUNS"])))
  assign[assign == 1] <- -1
  result <- data.frame(ptid = x$ptid[1], subset = unique(x$subset), assign = assign)
  return(result)
}

geomean <- function(x) exp(mean(log(x)))
library(flowReMix)
# Loading data -------------------------
tbdat <- readRDS("data/TB_rozot_booleans.rds")
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

configurations <- expand.grid(npost = c(10, 20, 40),
                              niter = c(20, 30, 40))
config <- configurations[setting, ]
npost <- as.numeric(config[1])
niter <- as.numeric(config[2])


# Analysis ----------------------------------
control <- flowReMix_control(updateLag = round(niter / 2), nsamp = 50, initMHcoef = 1,
                             nPosteriors = npost, centerCovariance = TRUE,
                             maxDispersion = 10^3, minDispersion = 10^7,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             lastSample = round(40 / npost), isingInit = -log(99),
                             initMethod = "robust",
                             preAssignCoefs = c(0.95, 0.5, seq(from = 0, to = 0.5, length.out = 10)),
                             ncores = ncores)

tempdat$stim <- tempdat$stimtemp
tempdat$stim[tempdat$stim == "UNS"] <- "aUNS"
tempdat$stim <- factor(tempdat$stim, levels = sort(unique(tempdat$stim)))
preAssign <- data.table::rbindlist(by(tempdat, tempdat$ptid, assign))
fit <- flowReMix(cbind(count, parentcount - count) ~ stim,
                 subject_id = ptid,
                 cell_type = subset,
                 cluster_variable = stim,
                 data = tempdat,
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "robust",
                 cluster_assignment = preAssign,
                 iterations = niter,
                 parallel = TRUE,
                 verbose = TRUE, control = control)

file <- paste("results/TBdat1_npost", npost, "_niter", niter, ".Robj", sep ="")
save(fit, file = file)





