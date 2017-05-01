#save(conditionalList, file = "data analysis/results/malariaConditional.Robj")
load(file = "data analysis/results/malariaConditional.Robj")

reps <- 50
resultlist <- list()
for(i in 1:length(conditionalList)) {
  mu <- conditionalList[[i]]$mle
  sigma <- conditionalList[[i]]$sigma
  threshold <- conditionalList[[i]]$threshold
  contrast <- conditionalList[[i]]$contrast
  result <- matrix(nrow = reps, ncol = 3)
  result <- data.frame(result)
  names(result) <- c("setting", "naive", "conditional")
  for(m in 1:reps) {
    pass <- FALSE
    while(!pass) {
      y <- as.numeric(mvtnorm::rmvnorm(1, mean = mu, sigma = sigma))
      product <- as.numeric(t(y) %*% contrast)
      if(product > threshold) pass <- TRUE
    }

    conditionalMLE <- linearConstraintMLE(y, sigma, contrast, threshold,
                                          alpha = 0.05,
                                          df = 1000,
                                          FWcorrection = TRUE)
    CI <- conditionalMLE$CI
    naive <- conditionalMLE$naiveCI
    cover <- numeric(length(y))
    naiveCover <- numeric(length(y))
    for(j in 1:nrow(CI)) {
      cover[j] <- CI[j, 1] < mu[j] & CI[j, 2] > mu[j]
      naiveCover[j] <- naive[j, 1] < mu[j] & naive[j, 2] > mu[j]
    }
    cover <- all(cover == 1)
    naiveCover <- all(naiveCover == 1)
    result[m, ] <- c(i, naiveCover, cover)
    print(names(conditionalList)[i])
    try(print(c(m, colMeans(result[1:m, 2:3, drop = TRUE]))))
  }

  resultlist[[i]] <- result
}
