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
cart <- readRDS("data/CART_stats_data3.rds", refhook = NULL)
names(cart) <- tolower(names(cart))
unique(cart$parent)
names(cart)[7] <- "date"
names(cart)[11] <- "isolation"
names(cart)[13] <- "dateTreated"
names(cart)[14] <- "plannedDose"
names(cart)[17] <- "crs"
names(cart)[22] <- "bestResponse"
cart <- subset(cart, !is.na(population))
by(cart, cart$id, function(x) unique(x$bestResponse))

cart <- cart[order(cart$id, cart$timepoint, cart$parent, cart$population), ]
cart[, c(5, 1, 4, 6, 22)]

# Collapse the IP timepoint
ip <- subset(cart, timepoint %in% c("CD4 IP", "CD8 IP"))
cart <- subset(cart, !(timepoint %in% c("CD4 IP", "CD8 IP")))
trimmed <- do.call("rbind", by(ip, ip$id, function(x) x[1, ]))
trimmed <- trimmed[, -c(1, 4, 2)]
trimmed$timepoint <- "IP"
ip <- summarize(group_by(ip, id, parent, population),
                  count = sum(count))
ip <- merge(ip, trimmed, all.x = TRUE)
cart <- rbind(cart, ip)

# Computing parentcount
parentcount <- summarize(group_by(cart, id, parent, timepoint),
                         parentcount = sum(count))
data.frame(parentcount)
cart <- merge(cart, parentcount, all.x = TRUE)
cart <- cart[order(cart$id, cart$timepoint, cart$parent, cart$population), ]
data.frame(cart$id, cart$population, cart$parentcount)

# Transforming timepoints to categories
cart$timenum <- 0
notip <- cart$timepoint != "IP"
times <- cart$timepoint[notip]
times <- as.numeric(substring(times, 2))
cart$timenum[notip] <- times
by(cart, cart$id, function(x) unique(x$timenum))
cart$timenum[cart$timenum > 0 & cart$timenum <= 19] <- 1
cart$timenum[cart$timenum > 19 & cart$timenum <= 35] <- 2
cart$timenum[cart$timenum > 35] <- 3
table(cart$timenum)

# Creating subset variable
cart$subset <- paste(cart$parent, cart$population, sep = "/")
by(cart, cart$id, function(x) unique(x$bestResponse))
by(cart, cart$subset, function(x) unique(x$poutcome))

# Exploration
unique(cart$population)
poutcome <- rep(TRUE, nrow(cart))
poutcome[cart$bestResponse %in% c("SD", "PD")] <- FALSE
poutcome[is.na(cart$bestResponse)] <- NA
cart$prop <- cart$count / cart$parentcount
cart$poutcome <- poutcome
# ggplot(cart) +
#   geom_boxplot(aes(x = factor(timenum), color = poutcome,
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


control <- flowReMix_control(updateLag = 5, nsamp = 50, initMHcoef = 1,
                             nPosteriors = 2, centerCovariance = TRUE,
                             maxDispersion = 10^8, minDispersion = 10^8,
                             randomAssignProb = 0.35, intSampSize = 50,
                             initMethod = "sparse")
system.time(fit <- flowReMix(cbind(count, parentcount - count) ~
                               treatment*timenum + disease + isolation + age + plannedDose,
                             cell_type = subset,
                             subject_id =  id, data = cart,
                             cluster_variable = treatment,
                             weights = NULL,
                             covariance = "sparse",
                             ising_model = "sparse",
                             regression_method = "sparse",
                             iterations = 10,
                             control = control))



