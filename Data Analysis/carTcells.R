library(flowReMix)
library(ggplot2)

# Processing Data -------------------
cart <- readRDS("data/merged_stats.rds", refhook = NULL)
names(cart) <- tolower(names(cart))
cart <- subset(cart, parent == "CAR-T")
by(cart, INDICES = cart$id, function(x) sort(unique(x$day)))
by(cart, INDICES = cart$id, function(x) unique(x$neurotox))
by(cart, INDICES = cart$id, function(x) unique(x$lymphodepletion))

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

subset(cart[, c(2, 13, 6, 16, 14)], id == "X329")

cart <- cart[order(cart$id, cart$dayfactor, cart$population), ]
remove <- by(cart, cart$population, function(x) length(unique(x$dayfactor)))
remove <- names(remove)[sapply(remove, function(x) x) == 1]
cart <- subset(cart, !(population %in% remove))

remove <- by(cart, INDICES = cart$id, function(x) sort(unique(x$day)))
remove <- names(remove)[which(sapply(remove, length) == 1)]
cart <- subset(cart, !(id %in% remove))

# Data exploration ---------------------------
# ggplot(cart) + geom_boxplot(aes(x = dayfactor, y = proportion, col = or_binary)) +
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
cartrep <- by(cart, cart$id, function(x) {
  x$weight <- 0.1
  x$treatment <- 0
  x$id <- runif(1)
  return(x)
  })
cartrep <- do.call("rbind", cartrep)
cart$treatment <- 1
cart$weight <- 1
cart <- rbind(cart, cartrep)
cart$population <- factor(cart$population)
cart$dayfactor <- as.character(as.numeric(cart$dayfactor))
system.time(fit <- subsetResponseMixtureRcpp(count ~  treatment*dayfactor,
                                             sub.population = population,
                                             N = parentcount, id =  id,
                                             data = cart,
                                             treatment = treatment,
                                             randomAssignProb = 0,
                                             weights = weight,
                                             updateLag = 5, nsamp = 50, maxIter = 15,
                                             initMHcoef = 3,
                                             covarianceMethod = "sparse",
                                             isingMethod = "sparse",
                                             regressionMethod = "sparse",
                                             centerCovariance = FALSE,
                                             dataReplicates = 5,
                                             maxDispersion = 10^3))



