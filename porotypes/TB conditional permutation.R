compMeasure <- function(TBlabels, posteriors) {
  signs <- colMeans(posteriors[TBlabels, ]) > colMeans(posteriors[-TBlabels, ])
  signs <- as.numeric(1 - 2 * (1 - signs))
  func <- posteriors %*% signs
  wilcoxTest <- wilcox.test(func[TBlabels], func[-TBlabels], exact = FALSE)$statistic
  return(wilcoxTest)
}

# aggregate <- readRDS(file = "data analysis/results/TBagg.rds")

posteriors <- aggregate$posteriors[, -1] %>% as.matrix()
tbdat <- aggregate$data
labels <- data.frame(ptid = tbdat$ptid, type = tbdat$type) %>% unique()
labels <- labels$type
nPerm <- 10^4
nTB <- sum(labels == "TB")
permutations <- replicate(nPerm, sample.int(length(labels), nTB, FALSE)) %>% t()
result <- apply(permutations, 1, compMeasure, posteriors)
observed <- compMeasure(which(labels == "TB"), posteriors)
hist(result)
abline(v = observed)
mean(observed <= result)
