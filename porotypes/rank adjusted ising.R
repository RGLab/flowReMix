# Setting up data from HVTN Fit --------------
load(file = "Data Analysis/results/HVTN bool betabinom 6 dropmore.Robj")
assignments <- fit$assignmentList
names(assignments) <- substr(names(assignments), 1, 9)
assignments <- lapply(unique(names(assignments)), function(x) {
  do.call("rbind", assignments[names(assignments) == x])
})
subsets <- names(fit$coefficients)
isCMV <- substring(subsets, 1, 3) == "cmv"
assignments <- lapply(assignments, function(x) x[, !isCMV])
subsets <- subsets[!isCMV]
expressed <- sapply(subsets, getExpression)
map <- cbind(subsets, expressed)
nsubsets <- ncol(assignments[[1]])
countCovar <- matrix(0, nrow = nsubsets, ncol = nsubsets)
mat <- t(sapply(assignments, function(x) x[sample(1:nrow(x), 1), ]))
colnames(mat) <- expressed
#ising <- IsingFit::IsingFit(mat, AND = TRUE, plot = TRUE)

# Fitting model --------------
nvars <- ncol(mat)
modelprobs <- (1 + nvars)^-1 / choose(nvars, 0:nvars)
offsets <- diff(log(modelprobs))
sum(modelprobs * choose(nvars, 0:nvars))
gamma <- 0.25
doParallel::registerDoParallel(cores=3)
minprob <- 1 / length(y)

system.time(isingmat <- foreach::foreach(j = 1:ncol(mat), .combine = rbind) %dopar% {
  print(j)
  y <- mat[, j]
  X <- mat[, -j]
  if(sum(y == 0) < 4) {
    p <- min(mean(y), 1 - minprob)
    row <- rep(0, ncol(mat))
    coef <- log(p / (1 - p))
    row[j] <- coef
    return(row)
  } else if(sum(y == 1) < 4) {
    p <- max(mean(y), minprob)
    row <- rep(0, ncol(mat))
    coef <- log(p / (1 - p))
    row[j] <- coef
    return(row)
  }
  off <- offsets[rowSums(X) + 1]
  netfit <- glmnet::glmnet(X, y, family = "binomial", offset = off)
  logliks <- 2 * (netfit$dev.ratio - 1) * netfit$nulldev
  dfs <- netfit$df
  ebic <- -logliks + dfs * log(nrow(mat) * (ncol(mat) - 1)^gamma)
  lambda <- netfit$lambda[which.min(ebic)]
  matrow <- rep(0, ncol(mat))
  coefs <- coef(netfit, s = lambda)
  matrow[-j] <- coefs[-1]
  matrow[j] <- coefs[1]
  return(matrow)
})
nonzero <- which(isingmat != 0, arr.ind = TRUE)
nonzero <- nonzero[which(nonzero[, 1] != nonzero[, 2]), ]
nonzero <- t(apply(nonzero, 1, sort))
nonzero <- unique(nonzero)
AND <- TRUE
for(i in 1:nrow(nonzero)) {
  u <- nonzero[i, 1]
  v <- nonzero[i, 2]
  first <- isingmat[u, v]
  second <- isingmat[v, u]
  if(AND & (first == 0 | second == 0)) {
    isingmat[u, v] <- 0
    isingmat[v, u] <- 0
    next
  }

  meanval <- (first + second) / 2
  isingmat[u, v] <- meanval
  isingmat[v, u] <- meanval
}



