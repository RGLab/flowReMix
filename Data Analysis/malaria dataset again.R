library(flowReMix)
# Malaria dataset ----------------------------
data(malaria)
names(malaria)
# table(malaria$experiment)
# unique(malaria$ptid)
# unique(malaria$population)
populations <- unique(malaria$population)
parents <- unique(malaria$parent)
leaves <- populations[!(populations %in% parents) ]
malaria <- subset(malaria, population %in% leaves)
# unique(malaria$stim)
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
control <- flowReMix_control(updateLag = 30, keepEach = 5, nsamp = 50, initMHcoef = 2,
                             nPosteriors = 4, centerCovariance = TRUE,
                             maxDispersion = 500, minDispersion = 10^8,
                             randomAssignProb = 10^-8, intSampSize = 100,
                             lastSample = 10, isingInit = -log(95),
                             initMethod = "robust",
                             ncores = 6,isingWprior=TRUE,
                             markovChainEM = TRUE,
                             prior = 0,
                             zeroPosteriorProbs = FALSE)

tempdat <- subset(malaria, TRUE)
tempdat$time <- tempdat$visitno
tempdat$subset <- factor(as.character(tempdat$subs))
tempdat$stim[tempdat$stim]
tempdat$trt <- 1
tempdat$visitno <- factor(as.character(tempdat$visitno), levels = c("Day 0", "Day 9", "pos", "Day 28", "Day 56", "Day 168"))
tempdat$visitInter <- tempdat$visitno
tempdat$trtTime <- tempdat$visitno
tempdat$stimInd <- as.numeric(tempdat$stim == "stim")

# tempdat <- subset(tempdat, parent == "4+" & stimgroup == "RBC")
# tempdat$subset <- factor(as.character(tempdat$subset))

# a <- model.matrix(cbind(count, parentcount - count) ~ visitno + stimInd:trtTime,
#              data = tempdat)
# colnames(a)

keep_levels = unlist(lapply(strsplit(levels(tempdat$subset),"/"),function(x)x[3]))%in%c("154+","IL21+","IL2+","IFNg+","TNFa+","IL4+")
temp = tempdat%>%filter(stimgroup=="RBC")%>%filter(subset %in% levels(subset)[keep_levels])%>%mutate(subset=factor(subset))

system.time(fit_model4_50 <- flowReMix(cbind(count, parentcount - count) ~  visitno + stimInd + stimInd:trtTime,
                 subject_id = ptid,
                 cell_type = subset,
                 cluster_variable = trtTime,
                 data = temp,
                 covariance = "sparse",
                 ising_model = "sparse",
                 regression_method = "robust",cluster_assignment = TRUE,
                 iterations = 100,
                 parallel = TRUE,
                 verbose = FALSE, control = control))

# load(file = "data analysis/results/malaria3_6_npost10niter40maxDisp1000.Robj")
# load(file = "data analysis/results/malaria3_5_npost5niter40maxDisp1000.Robj")
# load(file = "data analysis/results/malaria3_4_npost40niter30maxDisp1000.Robj")
# load(file = "data analysis/results/malaria3_3_npost20niter30maxDisp1000.Robj")
# load(file = "data analysis/results/malaria3_2_npost5niter20maxDisp100.Robj")
# load(file = "data analysis/results/malaria3_8_npost40niter40maxDisp1000.Robj")
# load(file = "data analysis/results/malaria3_7_npost20niter40maxDisp1000.Robj")

filenames <- as.list(dir(path = 'data analysis/results', pattern="malaria4visitno_*"))
filenames <- lapply(filenames, function(x) paste0('data analysis/results/', x))

filenames <- as.list(dir(path = 'data analysis/results', pattern="malaria3_*"))
filenames <- lapply(filenames, function(x) paste0('data analysis/results/', x))


add_ptid <- function(x, subject_id) {
  x$subject_id <- match.call()$subject_id
  return(x)
}

post <- list()
for(i in 1:length(filenames)) {
  load(file = filenames[[i]])
  post[[i]] <- fit$posteriors[, -1]
}
post <- Reduce("+", post) / length(filenames)
fit$posteriors[, -1] <- post
fit$data <- tempdat
fit <- add_ptid(fit, ptid)

# ROC table -------------------------------
outcome <- by(tempdat, tempdat$ptid, function(x) x$infection[1])
outcome <- data.frame(ptid = names(outcome), infection = as.numeric(outcome))
infection <- outcome[, 2]

rocResult <- rocTable(fit, target = infection, pvalue = "wilcoxon")
select <- apply(fit$posteriors[, -1], 2, function(x) max(x) >= 1 - 1)
select <- sapply(fit$coefficients, function(x) x[7] > 0)
select <- rep(TRUE, length(fit$coefficients))
rocResult$qvalue <- NA
rocResult$qvalue[select] <- p.adjust(rocResult$pvalue[select], method = "BH")
(rocResult[order(rocResult$auc, decreasing = TRUE), ])

# Plotting raw graphs ------------------
plot(fit, type = "graph", graph = "ising", threshold = 0.95,
     count = FALSE, fill = rocResult$auc)

# Boxplots ---------------------------
library(cowplot)
weights <- list()
weights$select <- select
weights$select <- rep(1, length(select))
boxall <- plot(fit, type = "boxplot", target = infection, groups = "all",
     test = "wilcoxon", one_sided = TRUE, jitter = TRUE)
boxall
# save_plot(boxall, filename = "figures/malariaBoxAll2.pdf",
#           base_height = 4, base_width = 5)

subsets <- colnames(fit$posteriors)[-1]
stims <- sapply(subsets, function(x) strsplit(x, "/")[[1]][1])
stimnames <- unique(stims)
stimgroups <- lapply(stimnames, function(x) subsets[stims == x])
names(stimgroups) <- stimnames
boxstim <- plot(fit, type = "boxplot", target = infection, test = "wilcoxon",
     groups = stimgroups, weights = weights, one_sided = TRUE, jitter = TRUE,
     ncol = 2)
boxstim
# save_plot(boxstim, filename = "figures/malariaBoxStim2.pdf",
#           base_height = 4, base_width = 7)


parents <- sapply(subsets, function(x) strsplit(x, "/")[[1]][2])
parentnames <- unique(parents)
parentgroups <- lapply(parentnames, function(x) subsets[parents == x])
names(parentgroups) <- parentnames
parentbox <- plot(fit, type = "boxplot", target = infection,
     groups = parentgroups,
     test = "wilcoxon", one_sided = TRUE, ncol = 3,
     jitter = TRUE)
parentbox
# save_plot(parentbox, filename = "figures/malariaBoxParent2.pdf",
#           base_height = 6, base_width = 7)

sc <- sapply(subsets, function(x) paste(strsplit(x, "/")[[1]][1:2], collapse = "/"))
scnames <- unique(sc)
scgroups <- lapply(scnames, function(x) subsets[sc == x])
names(scgroups) <- scnames
scbox <- plot(fit, type = "boxplot", target = infection,
     groups = scgroups,
     test = "wilcoxon", one_sided = TRUE, ncol = 4, jitter = TRUE)
scbox
# save_plot(scbox, filename = "figures/malariaSCbox2.pdf",
#           base_height = 6, base_width = 10)

# Graphical Models --------------------------
load(file = "data analysis/results/malariaIsing10.Robj")
plot(stability, threshold = 0)

load(file = "data analysis/results/malariaRand10.Robj")
plot(rand, threshold = 0.5, fill = rocResult$auc)


# Plotting coefficients and proportions -----------------
nplot <- 10
tempdat$prop <- tempdat$count / tempdat$parentcount
forplot <- tempdat[, c(9, 8, 10, 13, 14, 20)]
control <- subset(forplot, stim == "control")
control$stim <- NULL
names(control)[ncol(control)] <- "control"
stim <- subset(forplot, stim == "stim")
stim$stim <- NULL
names(stim)[ncol(stim)] <- "stim"
forplot <- merge(control, stim)
rocResult <- rocTable(fit, target = infection, pvalue = "wilcoxon",
                      sort = TRUE)
forplot <- subset(forplot, subset %in% rocResult$subset[1:nplot])
coefs <- fit$coefficients[names(fit$coefficients) %in% rocResult$subset[1:nplot]]
coefdat <- list()
expit <- function(x) 1 / (1 + exp(-x))
for(i in 1:length(coefs)) {
  coef <- coefs[[i]]
  dat <- data.frame(subset = names(coefs[i]),
                    visitno = c("Day 0", "Day 9", "pos", "Day 28", "Day 56", "Day 168"),
                    predcontrol = expit(sum(coef[1])),
                    predstim = c(expit(sum(coef[1:2])),
                                 expit(sum(coef[1:3])),
                                 expit(sum(coef[c(1:2, 4)])),
                                 expit(sum(coef[c(1:2, 5)])),
                                 expit(sum(coef[c(1:2, 6)])),
                                 expit(sum(coef[c(1:2, 7)]))))
  datnonresp <- data.frame(subset = names(coefs[i]),
                           visitno = c("Day 0", "Day 9", "pos", "Day 28", "Day 56", "Day 168"),
                           predcontrolnonresp = expit(sum(coef[1])),
                           predstimnonresp = c(expit(sum(coef[1:2])),
                                               expit(sum(coef[1:2])),
                                               expit(sum(coef[c(1:2)])),
                                               expit(sum(coef[c(1:2)])),
                                               expit(sum(coef[c(1:2)])),
                                               expit(sum(coef[c(1:2)]))))
  coefdat[[i]] <- merge(dat, datnonresp)
}
coefdat <- do.call("rbind", coefdat)
forplot <- merge(forplot, coefdat)
forplot <- subset(forplot, subset != "RBC/8+/IL21+")

size <- 0.4
ggplot(forplot, aes(x = as.numeric(visitno), y = stim - control, col = infection)) +
  geom_point(size = 0.5) + facet_wrap(~ subset, scale = "free_y") +
  geom_hline(yintercept = 0) + theme_bw() +
  geom_line(aes(y = predstim - predcontrol), col = "blue", size = size) +
  geom_line(aes(y = predstimnonresp - predcontrolnonresp), col = "red", linetype = 2, size = size)# +
  scale_colour_manual(values = c("red", "blue"))





