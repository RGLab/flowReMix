library(flowReMix)
library(ggplot2)
library(dplyr)
library(reshape2)
library(pROC)

activation <- c("CD69+", "CD95+", "CD122+", "HLA-DR+",
                "CD127+", "CD137+", "CD275+")
exhaustion <- c("CD160+", "CD244+", "CD25+", "KLRG1+",
                "LAG3+", "PD-1+", "TIGIT+", "TIM3+",
                "TIM3+")
tcell <- c("CD4+", "CD8+", "CD45RA+", "CD45RO+")
memory <- c("CCR7+", "RA-CCR7-", "RA-CCR7+", "RA+CCR7-",
            "RA+CCR7+", "RO+CCR7-", "RO+CCR7+", "CD27+",
            "CD28+", "CD45RA+", "CD45RO-", "CD95+")

# Processing Data -------------------
cart <- readRDS("data/carTstats.RDS", refhook = NULL)
names(cart) <- tolower(names(cart))
unique(cart$parent)
names(cart)[6] <- "ptid"
names(cart)[8] <- "date"
names(cart)[11] <- "isolation"
names(cart)[13] <- "dateTreated"
names(cart)[14] <- "plannedDose"
names(cart)[17] <- "crs"
names(cart)[18] <- "peakabsCD8"
names(cart)[19] <- "peakabsCD4"
names(cart)[20] <- "peakpercentCD8"
names(cart)[21] <- "peakpercentCD4"
names(cart)[22] <- "best"
cart <- subset(cart, !is.na(population))
timepoints <- by(cart, cart$ptid, function(x) unique(x$timepoint))
times <- unique(cart$timepoint)

cart <- cart[order(cart$ptid, cart$timepoint, cart$parent, cart$population), ]
cart[, c(5, 1, 4, 6, 22)]
cart$outcome <- "negative"
cart$outcome[cart$best %in% c("CR", "PR")] <- "positive"
cart$outcome[is.na(cart$best) | cart$best == "NA"] <- NA
cart <- subset(cart, cart$parentcount > 0)

# Keeping only leaves
leaves <- unique(cart$population)
parents <- unique(cart$parent)
leaves <- leaves[!(leaves %in% parents)]
cart <- subset(cart, population %in% leaves)

# Collapse the IP timepoint
ip <- subset(cart, timepoint %in% c("CD4 IP", "CD8 IP"))
cart <- subset(cart, !(timepoint %in% c("CD4 IP", "CD8 IP")))
trimmed <- do.call("rbind", by(ip, ip$ptid, function(x) x[1, ]))
trimmed <- trimmed[, -c(1, 4, 2)]
trimmed$timepoint <- "d0"
ip <- summarize(group_by(ip, ptid, parent, population),
                  count = sum(count), parentcount = sum(parentcount))
ip <- merge(ip, trimmed, all.x = TRUE)
ip$timepoint <- "d0"
cart <- rbind(cart[, -1], ip)

# Transforming timepoints to categories
cart$timenum <- 0
notip <- (cart$timepoint != "IP") | is.na(cart$timepoint)
times <- cart$timepoint[notip]
times <- as.numeric(substring(times, 2))
cart$timenum[notip] <- times
by(cart, cart$ptid, function(x) unique(x$timenum))
cart$timenum[cart$timenum > 0 & cart$timenum <= 19] <- 1
cart$timenum[cart$timenum > 19 & cart$timenum <= 35] <- 2
cart$timenum[cart$timenum > 35] <- 3
table(cart$timenum)

# Creating subset variable
cart <- subset(cart, !is.na(timenum))
cart$subset <- paste(cart$parent, cart$population, sep = "/")
#cart$subset <- paste(cart$timenum, cart$subset, sep = "/")
by(cart, cart$ptid, function(x) unique(x$best))
by(cart, cart$subset, function(x) unique(x$outcome))

# Exploration
unique(cart$population)
poutcome <- rep(TRUE, nrow(cart))
poutcome[cart$bestResponse %in% c("SD", "PD")] <- FALSE
poutcome[is.na(cart$bestResponse)] <- NA
cart$prop <- cart$count / cart$parentcount
cart$poutcome <- poutcome

# ggplot(cart) +
#   geom_boxplot(aes(x = factor(timenum), color = outcome,
#                    y = log(prop + 1 / parentcount))) +
#   facet_wrap(~ subset, scales = "free")

# Analysis
cart$treatment <- rep(1, nrow(cart))
cart$plannedDose[is.na(cart$plannedDose)] <- 4
cart$disease[is.na(cart$disease)] <- "missing"
cart$age[is.na(cart$age)] <- mean(cart$age[!is.na(cart$age)])
cart$isolation[is.na(cart$isolation)] <- "missing"
cart <- subset(cart, population != "!LAG3&!PD1&!TIM3")
cart$subset <- factor(cart$subset)
cart <- subset(cart, timenum < 4)
cart$timenum <- factor(cart$timenum)


library(flowReMix)
control <- flowReMix_control(updateLag = 3, nsamp = 100, initMHcoef = 1,
                             nPosteriors = 1, centerCovariance = TRUE,
                             maxDispersion = 10^2 * 5, minDispersion = 10^4,
                             randomAssignProb = 10^-8, intSampSize = 50,
                             lastSample = 40, isingInit = -log(89),
                             initMethod = "robust")

tempdat <- subset(cart, !is.na(outcome))
system.time(fit <- flowReMix(cbind(count, parentcount - count) ~ timenum*outcome,
                             subject_id = ptid,
                             cell_type = subset,
                             cluster_variable = NULL,
                             data = tempdat,
                             covariance = "sparse",
                             ising_model = "sparse",
                             regression_method = "robust",
                             iterations = 6,
                             parallel = TRUE,
                             verbose = TRUE, control = control))
#save(fit, file = "data analysis/results/carTagain 3.Robj")
#load(file = "data analysis/results/carTagain 1.Robj")
post <- fit$posteriors
post <- merge(post, unique(data.frame(ptid = cart$ptid, outcome =  cart$outcome)),
              all.x = TRUE, all.y = FALSE)
post <- post[, -1]
outcome <- post[, ncol(post)]
post <- post[, -ncol(post)]
#outcome[is.na(outcome)] <- "negative"
library(pROC)
aucs <- numeric(ncol(post))
for(i in 1:ncol(post)) {
  aucs[i] <- roc(outcome ~ post[, i])$auc
}
n0 <- sum(outcome == "negative", na.rm = TRUE)
n1 <- sum(outcome == "positive", na.rm = TRUE)
pvals <- pwilcox(aucs * n0 * n1, n0, n1, lower.tail = FALSE)
results <- data.frame(probs = fit$levelProbs, auc = aucs, pval = pvals)
results <- results[order(results$auc), ]
results <- subset(results, probs > 0.05 & probs < 0.95)
qvals <- p.adjust(results$pval, method = "BH")


# Graphical Model -------------
graph <- fit$isingfit
diag(graph) <- 0
which(graph != 0, arr.ind = TRUE)
library(igraph)
graph <- graph.adjacency(graph)
comp <- components(graph)
aucs <- numeric(2)
g <- 1
test <- which(comp$csize > 4)
for(i in c(test)) {
  group <- which(comp$membership == i)
  score <- rowMeans(post[, group])
  aucs[g] <- roc(outcome ~ score)$auc
  g <- g + 1
}
n0 <- sum(outcome == "negative", na.rm = TRUE)
n1 <- sum(outcome == "positive", na.rm = TRUE)
pvals <- pwilcox(aucs * n0 * n1, n0, n1, lower.tail = FALSE)

