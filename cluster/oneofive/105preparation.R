filtered_data <- readRDS(file="data/hvtn105_filtered_data.rds")

# Defining (time X stim X cell) subsets
names(filtered_data) <- tolower(names(filtered_data))
unique(filtered_data$stim)
filtered_data <- filtered_data[, -14]
filtered_data$subset <- gsub(".", "", filtered_data$subset, fixed = "TRUE")
stimdat <- stimulationModel(filtered_data,
                            cell_type = subset,
                            stim_var = stim,
                            controls = "negctrl",
                            stim_groups = list(gag = "Gag-ZM96",
                                               vacenv = "92TH023-ENV",
                                               env = c("Env-1-ZM96", "Env-2-ZM96")))
unique(stimdat$stimCellType)
stimdat$visitStimCell <- interaction(stimdat$visitno, stimdat$stimCellType, sep = "/", drop = TRUE)

# Screening samples with very few cells ----
stimdat$totalCells <- stimdat$name %>% strsplit(split = "_") %>%
  sapply(function(x) x[2]) %>% as.numeric()
stimdat <- subset(stimdat, totalCells > 10^5)

# Screening samples with very small parent count ----
stimdat <- subset(stimdat, parentcount >= 100)

# Screening samples that are not paired ------
paired <- group_by(stimdat, visitStimCell, ptid) %>%
  summarize(pairedInd = sum(stim == "ctrl") > 0 & sum(stim != "ctrl") > 0)
stimdat <- merge(stimdat, paired, all.x = TRUE)
stimdat <- subset(stimdat, pairedInd)

# Screening based on counts:
# Keeping cell-subsets that have at least 10 subjects with non-zero counts in all time points
countSummary <- group_by(stimdat, stimCellType, visitno, ptid) %>%
  summarize(nonzero = any(count > 0)) %>%
  group_by(stimCellType, visitno) %>% summarize(nonzero = sum(nonzero)) %>%
  group_by(stimCellType) %>% summarize(nonzero = mean(nonzero >= 10)) %>%
  data.frame() %>% subset(nonzero == 1)

stimdat <- subset(stimdat, stimCellType %in% countSummary$stimCellType)
stimdat$visitStimCell <- factor(as.character(stimdat$visitStimCell))
stimdat$stimCellType <- factor(stimdat$stimCellType)
levels(stimdat$visitStimCell)
stimdat <- data.frame(stimdat)

# Number of subjects with each subset ------------
subsets <- stimdat$stimCellType %>% levels()
subjectsPerSub <- group_by(stimdat, stimCellType, visitno) %>%
  summarize(nSubjects = length(unique(ptid))) %>%
  data.table::setorder(nSubjects) %>%
  data.frame() %>% dcast( stimCellType ~ visitno  , value.var = "nSubjects")
names(subjectsPerSub)[1] <- "subset"
subjectsPerSub$subset <- as.character(subjectsPerSub$subset)

# Screening with lme4 -------
# library(lme4)
# subsets <- stimdat$stimCellType %>% levels()
# screenResults <- data.frame(subset = subsets, pval = 1)
# for(i in 1:length(subsets)) {
#   pop <- screenResults$subset[i] %>% as.character()
#   subdat <- subset(stimdat, stimCellType == pop)
#   subdat$od <- 1:nrow(subdat)
#   lmeFit <- NULL
#   try(lmeFit <- glmer(cbind(count, parentcount - count) ~ stim*visitno + (1|od) + (1|ptid),
#                   data = subdat, family = "binomial") %>% summary())
#   if(is.null(lmeFit)) {
#     next
#   }
#   nullFit <- NULL
#   try(nullFit <- glmer(cbind(count, parentcount - count) ~ visitno + (1|od) + (1|ptid),
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
# saveRDS(screenResults, file = "data/HVTN105screen.rds")
qvalThreshold <- 0.05
targetSize <- 110
minGroupSize <- 4
minSubjectCount <- 89

screenResults <- readRDS(file = "data/HVTN105screen.rds")
screenResults$subset <- as.character(screenResults$subset)
screenResults <- merge(screenResults, subjectsPerSub)
screenResults <- screenResults[order(screenResults$pval), ]
screenResults$minCount <- apply(screenResults[, 3:5], 1, min)
screenResults$minCount[is.na(screenResults$minCount)] <- 0
screenResults <- subset(screenResults, minCount >= minSubjectCount)
screenResults$qval <- p.adjust(screenResults$pval, method = "BH")

sapply(c(10^-6, 10^-5, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.15, 0.2), function(x) sum(screenResults$qval <= x))
screenResults[1:102, ]
screenResults$stimpar <- screenResults$subset %>% as.character() %>%
  strsplit("/", fixed = TRUE) %>%
  sapply(function(x) paste(x[1:2], collapse = "/"))
groupScreen <- by(screenResults, screenResults$stimpar, function(x) sum(x$qval < qvalThreshold)) %>%
  c()
groupScreen <- groupScreen[groupScreen >= minGroupSize]
keepGroups <- names(groupScreen)
screenResults <- subset(screenResults, qval < qvalThreshold & stimpar %in% keepGroups)
groupedScreen <- split(screenResults, screenResults$stimpar)
selectedSubsets <- vector(nrow(screenResults), mode = "character")
nextInd <- 1
topN <- 1
while(nextInd <= targetSize) {
  for(i in 1:length(groupedScreen)) {
    if(nrow(groupedScreen[[i]]) >= topN) {
      selectedSubsets[nextInd] <- groupedScreen[[i]]$subset[topN]
      nextInd <- nextInd + 1
    }
  }
  topN <- topN + 1
}
selectedSubsets <- selectedSubsets[!(selectedSubsets == "")]
screenResults <- subset(screenResults, subset %in% selectedSubsets)
table(screenResults$stimpar)
nrow(screenResults)
stimdat <- subset(stimdat, stimCellType %in% selectedSubsets)
# saveRDS(stimdat, file = "data/hvtn105_110subsets.rds")




