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
control <- flowReMix_control(nsamp = 15, dataReplicates = 3,
                             centerCovariance = FALSE, updateLag = 5)
system.time(fit <- flowReMix(cbind(count, parentcount - count) ~ treatment + icu + disease,
                             cell_type = subset,
                             subject_id =  id, data = cart,
                             cluster_variable = treatment,
                             weights = NULL,
                             covariance = "sparse",
                             ising_model = "sparse",
                             regression_method = "sparse",
                             iterations = 10,
                             control = control))

#save(fit, file = "data analysis/results/carT results.Robj")
require(pROC)
par(mfrow = c(3, 5), mar = rep(1, 4))
subsets <- names(fit$posteriors)[-1]
pvals <- numeric(length(subsets))
aucs <- pvals
for(i in 1:length(subsets)) {
  rocfit <- roc(outcome ~ fit$posteriors[, i + 1])
  auc <- rocfit$auc
  n1 <- sum(outcome == 1)
  n2 <- sum(outcome == 2)
  Ustat <- auc * n1 * n2
  pval <- min(1, 2 * (1 - pwilcox(Ustat, n1, n2)))
  pvals[i] <- pval
  aucs[i] <- auc
  plot(rocfit, main = paste(paste(subsets[i], " ", round(auc, 3), " ", round(pval, 3))))
}
dataks <- ks.test(pvals, punif, alternative = "greater")[[2]]

qvals <- p.adjust(pvals, method = "BH")
names(qvals) <- subsets
categories <- list(activation = pvals[subsets %in% activation],
          memory = pvals[subsets %in% memory],
          tcells = pvals[subsets %in% tcell],
          exhaustion = pvals[subsets %in% exhaustion])

categories <- lapply(categories, function(x) p.adjust(x, method = "BH"))

# wilcoxon distribution
p <- length(qvals)
wilcsamp <- rwilcox(10^4 * p, n1, n2)
wilcpvals <- pwilcox(wilcsamp, n1, n2)
wilcsamp <- 2 * pmin(wilcpvals, 1 - wilcpvals)
wilcsamp <- matrix(wilcsamp, ncol = p)
wilcsamp <- t(apply(wilcsamp, 1, sort))
meanwilc <- colMeans(wilcsamp)
diff <- apply(wilcsamp, 1, function(x) mean((x - meanwilc)^2))
actualdiff <- mean((sort(pvals) - meanwilc)^2)
par(mfrow = c(1, 1))
hist(diff)
abline(v = actualdiff)
1 - mean(diff < actualdiff)

all <- data.frame(pvals, aucs, subsets)[order(pvals), ]
all$meanwilc <- meanwilc
all


# Better estimation of graphical model ----------------------
assignments <- fit$assignmentList
names(assignments) <- sapply(names(assignments), stringr::str_sub, 1, 4)
assignments <- lapply(unique(names(assignments)),
                      function(x) do.call("rbind", assignments[names(assignments) == x]))
reps <- 1000
modelList <- list()
nsubsets <- ncol(assignments[[1]])
countCovar <- matrix(0, nrow = nsubsets, ncol = nsubsets)
for(i in 1:reps) {
  mat <- t(sapply(assignments, function(x) x[sample(1:nrow(x), 1), ]))
  keep <- apply(mat, 2, function(x) any(x != x[1]))
  mat <- mat[, keep]
  model <- IsingFit::IsingFit(mat, AND = FALSE, plot = FALSE)
  modelList[[i]] <- model
  #plot(model)
  countCovar[keep, keep] <- countCovar[keep, keep] + (model$weiadj != 0) * sign(model$weiadj)
  print(i)
}

props <- countCovar / reps
table(props)
threshold <- 0.25
which(props > threshold, arr.ind = TRUE)
props[abs(props) <= threshold] <- 0
sum(props != 0) / 2

# Plotting graph ---------------------
require(GGally)
library(network)
library(sna)
network <- props
keep <- apply(network, 1, function(x) any(abs(x) > threshold))
network <- network[keep, keep]
net <- network(props)
nodes <- ggnet2(network, label = subsets[keep])$data
edges <- matrix(nrow = sum(network != 0)/2, ncol = 5)
p <- nrow(network)
row <- 1
for(j in 2:p) {
  for(i in 1:(j-1)) {
    if(network[i, j] != 0) {
      edges[row, ] <- unlist(c(nodes[i, 6:7], nodes[j, 6:7], network[i, j]))
      row <- row + 1
    }
  }
}

edges <- data.frame(edges)
names(edges) <- c("xstart", "ystart", "xend", "yend", "width")
coefficients <- sapply(fit$coefficients, function(x) x[2])
nodes$probs <- fit$levelProbs[keep]
nodes$auc <- aucs[keep]

names(edges)[5] <- "Dependence"
lims <- max(abs(props))
ggplot() +
  scale_colour_gradient2(limits=c(-lims, lims), low="dark red", high = "dark green") +
  geom_segment(data = edges, aes(x = xstart, y = ystart,
                                 xend = xend, yend = yend,
                                 col = Dependence,
                                 alpha = Dependence),
               size = 1) +
  scale_fill_gradient2(low = "white", high = "red", limits = c(0.5, 1)) +
  #scale_fill_gradientn(colours = rainbow(4))+
  geom_point(data = nodes, aes(x = x, y = y, fill = auc), shape = 21,
             size = 12, col = "grey") +
  scale_shape(solid = FALSE) +
  geom_text(data = nodes, aes(x = x, y = y, label = nodes$label), size = 1.8) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(), legend.position = "none")

nodes[order(nodes$auc),]


# Permutation test for permutation test!
m <- 400
p <- length(subsets)
n <- length(outcome)
kspvals <- numeric(m)
n1 <- sum(outcome == 1)
n2 <- sum(outcome == 2)
for(j in 1:m) {
  cat(j, ' ')
  aucs <- numeric(p)
  pvals <- numeric(p)
  tempoutcome <- outcome[order(runif(n))]
  for(i in 1:p) {
    rocfit <- roc(tempoutcome ~ fit$posteriors[, i + 1])
    auc <- rocfit$auc
    Ustat <- auc * n1 * n2
    pval <- min(1, 2 * (1 - pwilcox(Ustat, n1, n2)))
    pvals[i] <- pval
    aucs[i] <- auc
  }

  kspvals[j] <- ks.test(pvals, punif, alternative = "greater")[[2]]
}
mean(kspvals < dataks)


