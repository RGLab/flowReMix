library(flowReMix)
library(ggplot2)

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
cart <- readRDS("data/merged_stats.rds", refhook = NULL)
names(cart) <- tolower(names(cart))
cart <- subset(cart, parent == "CAR-T")
by(cart, INDICES = cart$id, function(x) sort(unique(x$day)))
by(cart, INDICES = cart$id, function(x) unique(x$neurotox))
by(cart, INDICES = cart$id, function(x) unique(x$lymphodepletion))
by(cart, INDICES = cart$id, function(x) c(unique(x$or_binary), length(unique(x$day))))

# standardizing population names
populations <- strsplit(cart$population, split = "/")
populations <- sapply(populations, function(x) x[[length(x)]])
cart$population <- populations
# setting days into categories
cart$dayfactor <- 0
cart$dayfactor[cart$day > 0 & cart$day <= 9] <- 1
cart$dayfactor[cart$day > 9 & cart$day <= 15] <- 2
cart$dayfactor[cart$day > 15 & cart$day <= 25] <- 3
cart$dayfactor[cart$day > 25] <- 4
cart$dayfactor <- factor(cart$dayfactor, levels = c(0:4))
cart$proportion <- cart$count / cart$parentcount
tempcart <- subset(cart, dayfactor == "0")
by(tempcart, tempcart$population, function(x) data.frame(x$id, x$count, x$proportion,
                                                         x$neurotox, x$lymphodepletion,
                                                         x$or_binary))

subset(cart[, c(2, 13, 6, 16, 14)], id == "X329")

# cart <- cart[order(cart$id, cart$dayfactor, cart$population), ]
# remove <- by(cart, cart$population, function(x) length(unique(x$dayfactor)))
# remove <- names(remove)[sapply(remove, function(x) x) == 1]
# cart <- subset(cart, !(population %in% remove))

# remove <- by(cart, INDICES = cart$id, function(x) sort(unique(x$day)))
# remove <- names(remove)[which(sapply(remove, length) == 1)]
# cart <- subset(cart, !(id %in% remove))

# Data exploration ---------------------------
cart$subset <- paste(cart$tcellpop, cart$population, sep = "/")
ggplot(subset(cart, dayfactor == "0"), aes(x = subset, y = log(proportion), col = or_binary)) +
         geom_boxplot() + geom_jitter(size = 0.1)
#
# ggplot(cart) +
#   geom_line(aes(x = day, y = proportion, col = or_binary, linetype = id)) +
#   facet_wrap(~ population, ncol = 6, scales = "free_y")

correlations <- by(cart, list(cart$id, cart$dayfactor), function(x) cbind(x$proportion))
correlations <- correlations[sapply(correlations, length) == 27]
correlations <- t(do.call("cbind", correlations))
correlations <- cor(correlations)
diag(correlations) <- 0
populations <- unique(cart$population)
highcorr <- which(correlations > 0.85, arr.ind = TRUE)
highcorr <- t(apply(highcorr, 1, sort))
highcorr <- unique(highcorr)
cbind(populations[highcorr[, 1]], populations[highcorr[, 2]])

# Screening cell-subsets?

# Attempt to fit a model
outcome <- as.vector(by(cart, cart$id, function(x) x$or_binary[1]))
neurotox <- as.vector(by(cart, cart$id, function(x) x$neurotox[1]))
lymph <- as.vector(by(cart, cart$id, function(x) x$lymphodepletion[1]))

#cart <- subset(cart, dayfactor == "0")
cart$population <- factor(cart$population)
cart$treatment <- rep(1, nrow(cart))
cart <- subset(cart, dayfactor == "0")
cart$subset <- as.factor(as.character(paste(cart$tcellpop, cart$population)))
system.time(fit <- subsetResponseMixtureRcpp(count ~  treatment*icu*tcellpop*disease + age,
                                             sub.population = population,
                                             N = parentcount, id =  id,
                                             data = cart,
                                             treatment = treatment,
                                             randomAssignProb = 0.0,
                                             weights = NULL,
                                             updateLag = 5, nsamp = 10,
                                             maxIter = 15,
                                             initMHcoef = 0.5,
                                             covarianceMethod = "sparse",
                                             isingMethod = "sparse",
                                             regressionMethod = "sparse",
                                             centerCovariance = FALSE,
                                             dataReplicates = 20,
                                             maxDispersion = 10^3))

require(pROC)
par(mfrow = c(3, 5), mar = rep(1, 4))
subsets <- names(fit$posteriors)[-1]
pvals <- numeric(length(subsets))
for(i in 1:length(subsets)) {
  rocfit <- roc(outcome ~ fit$posteriors[, i + 1])
  auc <- rocfit$auc
  n1 <- sum(outcome == 1)
  n2 <- sum(outcome == 2)
  Ustat <- auc * n1 * n2
  pval <- min(1, 2 * (1 - pwilcox(Ustat, n1, n2)))
  pvals[i] <- pval
  plot(rocfit, main = paste(paste(subsets[i], " ", round(auc, 3), " ", round(pval, 3))))
}

qvals <- p.adjust(pvals, method = "BH")
names(qvals) <- subsets
categories <- list(activation = pvals[subsets %in% activation],
          memory = pvals[subsets %in% memory],
          tcells = pvals[subsets %in% tcell],
          exhaustion = pvals[subsets %in% exhaustion])

categories <- lapply(categories, function(x) p.adjust(x, method = "BH"))

