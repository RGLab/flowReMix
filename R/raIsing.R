raIsing <- function(mat, AND = TRUE, gamma = 0.9,
                    modelprobs = NULL, minprob = NULL,
                    method = "sparse") {
  nvars <- ncol(mat)
  if(!is.null(modelprobs) & length(modelprobs) != (ncol(mat) + 1)) {
    warning("modelprobs must be of length ncol(mat) + 1 !")
  }
  if(is.null(modelprobs)) {
    modelprobs <- (1 + nvars)^-1 / choose(nvars, 0:nvars)
  }
  offsets <- diff(log(modelprobs))

  if(gamma < 0) gamma <- 0
  if(foreach::getDoParWorkers() == 1) {
    foreach::registerDoSEQ()
  }

  if(is.null(minprob)) {
    minprob <- 1 / nrow(mat)
  }

  isingmat <- foreach(j = 1:ncol(mat), .combine = rbind) %dopar% {
    print(j)
    y <- as.vector(mat[, j])
    X <- as.matrix(mat[, -j])
    xcols <- colSums(X)
    regX <- X[, xcols >= 4]

    if(sum(y == 0) < 8) {
      p <- min(mean(y), 1 - minprob)
      row <- rep(0, ncol(mat))
      coef <- log(p / (1 - p))
      row[j] <- coef
      return(row)
    } else if(sum(y == 1) < 8) {
      p <- max(mean(y), minprob)
      row <- rep(0, ncol(mat))
      coef <- log(p / (1 - p))
      row[j] <- coef
      return(row)
    } else if(ncol(X) < 2) {
      p <- mean(y)
      row <- rep(0, ncol(mat))
      log(p / (1 - p))
      row[j] <- coef
      return(row)
    }

    # if(method == "joint") {
    #   off <- rep(-3, length(y))
    # } else {
    # }
    off <- offsets[rowSums(X) + 1]
    netfit <- glmnet::glmnet(regX, y, family = "binomial", offset = off,
                             intercept = TRUE)
    logliks <- 2 * (netfit$dev.ratio - 1) * netfit$nulldev
    dfs <- netfit$df
    ebic <- -logliks + dfs * log(nrow(mat) * (ncol(mat) - 1)^gamma)
    lambda <- netfit$lambda[which.min(ebic)]
    matrow <- rep(0, ncol(mat))
    coefs <- rep(0, ncol(mat))
    coefs[c(1, which(xcols >= 4) + 1)] <- coef(netfit, s = lambda)
    matrow[-j] <- coefs[-1]
    matrow[j] <- coefs[1]
    return(matrow)
  }

  nonzero <- which(isingmat != 0, arr.ind = TRUE)
  nonzero <- nonzero[which(nonzero[, 1] != nonzero[, 2]), ]
  nonzero <- t(apply(nonzero, 1, sort))
  nonzero <- unique(nonzero)
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

  # if(method == "joint") {
  #   diag(isingmat) <- -3
  # }

  return(isingmat)
}
