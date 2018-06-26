logit <- function(x) log(x / (1 - x))
expit <- function(x) 1 / (1 + exp(-x))

raIsing <- function(mat, AND = TRUE, gamma = 0.9,
                    modelprobs = NULL, minprob = NULL,
                    method = "sparse", cv = FALSE,
                    family = "binomial",verbose=FALSE,
                    weights = NULL, parallel = FALSE,
                    returnPath = FALSE) {
  if(is.null(weights)) {
    weights <- rep(1, nrow(mat))
  }

  nvars <- ncol(mat)
  if(!is.null(modelprobs) & length(modelprobs) != (ncol(mat) + 1)) {
    warning("modelprobs must be of length ncol(mat) + 1 !")
  }
  if(is.null(modelprobs)) {
    modelprobs <- (1 + nvars)^-1 / choose(nvars, 0:nvars)
  }
  offsets <- diff(log(modelprobs))

  if(gamma < 0) gamma <- 0

  if(is.null(minprob)) {
    minprob <- 1 / nrow(mat)
  }

  # Estimating the ising parameters via neighborhood selection
  isingmat <- foreach(j = 1:ncol(mat)) %dorng% {
    getNeighborhood(j, mat, family, off, gamma, weights, cv, method, minprob,
                    returnPath = returnPath)
  }
  if(returnPath) {
    isingmat <- lapply(isingmat, function(x) {
      if(is.list(x)){
        return(x)
      } else {
        l <- list()
        l$matrow <- x
        l$orderRow <- rep(length(x), length(x))
        return(l)
      }
      })
    pathmat <- lapply(isingmat, function(x) x$orderRow) %>% do.call("rbind", .)
    isingmat <- lapply(isingmat, function(x) x$matrow) %>% do.call("rbind", .)
  } else {
    isingmat <- do.call("rbind", isingmat)
  }

  nonzero <- which(isingmat != 0, arr.ind = TRUE)
  nonzero <- nonzero[which(nonzero[, 1] != nonzero[, 2]), , drop = FALSE]
  if(length(nonzero) != 0) {
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
  }

  if(returnPath) {
    pathmat <- (pathmat + t(pathmat)) / 2
    return(list(isingmat = isingmat, pathmat = pathmat))
  } else {
    return(isingmat)
  }
}

pIsing <- function(mat, AND = TRUE, gamma = 0.9,
                   method = "sparse", cv = FALSE,
                   empBayes = FALSE, preAssignment,
                   family = "binomial", prevfit, verbose=FALSE,
                   weights = NULL) {
  if(is.null(weights)) {
    weights <- rep(1, nrow(mat))
  }
  nvars <- ncol(mat)

  if(gamma < 0) gamma <- 0

  # Computing Diagonal -----------
  minProb <- 1 / nrow(mat)
  isingOffset <- numeric(ncol(mat))
  targetResp <- numeric(ncol(mat))
  for(i in 1:ncol(mat)) {
    pResp <- mean(mat[, i])
    pResp <- min(pResp, 1 - mean(preAssignment[preAssignment[, 2] == colnames(mat)[i], 3] == 0))
    pResp <- max(pResp, mean(preAssignment[preAssignment[, 2] == colnames(mat)[i], 3] == 1))
    pResp <- pmin(pmax(minProb, pResp), 1 - minProb)
    targetResp[i] <- pResp

    expit <- function(x) 1 / (1 + exp(-x))
    coefs <- as.numeric(prevfit[i, -i])
    covs <- as.matrix(mat[, -i])
    eta <- as.numeric(covs %*% coefs)
    # if(verbose) cat(pResp, " ")
    isingOffset[i] <- uniroot(f = function(off) weightedMean(expit(eta + off), w = weights,na_rm=FALSE) - pResp,
                              interval = c(-50, 50))$root
  }
  # if(verbose) cat("\n")

  isingmat <- foreach(j = 1:ncol(mat), .combine = rbind) %dorng% {
    y <- as.vector(mat[, j])
    X <- as.matrix(mat[, -j])
    xcols <- colSums(X)
    if(family == "binomial") {
      regX <- X[, xcols >= 4,drop=FALSE]
    } else {
      regX <- X
    }

    row <- rep(0, ncol(mat))
    if(sum(y == 0) < 8 & family == "binomial") {
      row[j] <- isingOffset[j]
      return(row)
    } else if(sum(y == 1) < 8 & family == "binomial") {
      row[j] <- isingOffset[j]
      return(row)
    } else if(ncol(regX) < 2) {
      row[j] <- isingOffset[j]
      return(row)
    }

    off <- rep(isingOffset[j], nrow(regX))

    if(!cv) {
      netfit <- glmnet::glmnet(regX, y, family = family, offset = off,
                               intercept = FALSE, weights = weights)
      logliks <- 2 * (netfit$dev.ratio - 1) * netfit$nulldev
      dfs <- netfit$df
      ebic <- -logliks + dfs * log(nrow(mat) * (ncol(mat) - 1)^gamma)
      lambda <- netfit$lambda[which.min(ebic)]
    } else {
      netfit <- glmnet::cv.glmnet(regX, y, family = family, offset = off,
                                  intercept = FALSE, weights = weights)
      lambda <- netfit$lambda.min
    }
    matrow <- rep(0, ncol(mat))
    coefs <- rep(0, ncol(mat))
    if(family == "binomial") {
      coefs[c(1, which(xcols >= 4) + 1)] <- coef(netfit, s = lambda)
    } else {
      coefs <- coef(netfit, s = lambda)
    }
    matrow[-j] <- coefs[-1]
    matrow[j] <- isingOffset[j]
    return(matrow)
  }

  nonzero <- which(isingmat != 0, arr.ind = TRUE)
  nonzero <- nonzero[which(nonzero[, 1] != nonzero[, 2]), , drop=FALSE ]
  if(length(nonzero) != 0) {
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

    for(i in 1:ncol(isingmat)) {
      coefs <- as.numeric(isingmat[i, -i])
      covs <- as.matrix(mat[, -i])
      eta <- as.numeric(covs %*% coefs)
      isingOffset[i] <- uniroot(f = function(off) weightedMean(expit(eta + off), w = weights,na_rm=FALSE) - targetResp[i],
                                interval = c(-50, 50))$root
    }
    diag(isingmat) <- isingOffset
  }

  return(isingmat)
}

getNeighborhood <- function(j, mat, family, off, gamma, weights, cv, method, minprob,
                            returnPath = FALSE) {
  y <- as.vector(mat[, j])
  X <- as.matrix(mat[, -j,drop=FALSE])
  xcols <- colSums(X)
  if(family == "binomial") {
    regX <- as.matrix(X[, xcols >= 4,drop=FALSE])
  } else {
    regX <- X
  }

  if(sum(y == 0) < 8 & family == "binomial"| ncol(regX) < 2) {
    p <- min(mean(y), 1 - minprob)
    row <- rep(0, ncol(mat))
    coef <- log(p / (1 - p))
    row[j] <- coef
    return(row)
  } else if(sum(y == 1) < 8 & family == "binomial") {
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

  if(method == "raIsing") {
    off <- offsets[rowSums(X) + 1]
  } else {
    off <- NULL
  }

  if(!cv) {
    netfit <- glmnet::glmnet(regX, y, family = family, offset = off,
                             intercept = TRUE, weights = weights)
    logliks <- 2 * (netfit$dev.ratio - 1) * netfit$nulldev
    dfs <- netfit$df
    ebic <- -logliks + dfs * log(nrow(mat) * (ncol(mat) - 1)^gamma)
    lambda <- netfit$lambda[which.min(ebic)]
  } else {
    netfit <- glmnet::cv.glmnet(regX, y, family = family, offset = off,
                                intercept = TRUE, weights = weights)
    lambda <- netfit$lambda.min
  }
  matrow <- rep(0, ncol(mat))
  coefs <- rep(0, ncol(mat))
  if(family == "binomial") {
    coefs[c(1, which(xcols >= 4) + 1)] <- coef(netfit, s = lambda)
  } else {
    coefs <- coef(netfit, s = lambda)
  }
  matrow[-j] <- coefs[-1]
  matrow[j] <- coefs[1]

  if(returnPath) {
    optimPath <- apply(as.matrix(netfit$beta), 1, function(x) {
      if(any(x != 0)) {
        return(min(which(x != 0)))
      } else {
        return(length(matrow))
      }}) %>% order()
    orderRow <- rep(length(matrow), length(matrow))
    orderRow[j] <- 0
    orderRow[-j][which(xcols >= 4)][order(optimPath)] <- 1:sum(xcols >= 4)
    return(list(matrow = matrow, orderRow = orderRow))
  } else {
    return(matrow)
  }
}
