logit <- function(x) log(x / (1 - x))
expit <- function(x) 1 / (1 + exp(-x))

#' @export
raIsing <- function(mat, AND = TRUE, gamma = 0.9,
                    modelprobs = NULL, minprob = NULL,
                    method = "sparse", cv = FALSE,
                    family = "binomial",verbose=FALSE) {
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
    y <- as.vector(mat[, j])
    X <- as.matrix(mat[, -j])
    xcols <- colSums(X)
    if(family == "binomial") {
      regX <- X[, xcols >= 4]
    } else {
      regX <- X
    }

    if(sum(y == 0) < 8 & family == "binomial") {
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
                               intercept = TRUE)
      logliks <- 2 * (netfit$dev.ratio - 1) * netfit$nulldev
      dfs <- netfit$df
      ebic <- -logliks + dfs * log(nrow(mat) * (ncol(mat) - 1)^gamma)
      lambda <- netfit$lambda[which.min(ebic)]
    } else {
      netfit <- glmnet::cv.glmnet(regX, y, family = family, offset = off,
                                  intercept = TRUE)
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
    return(matrow)
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

  return(isingmat)
}

pIsing <- function(mat, AND = TRUE, gamma = 0.9,
                   method = "sparse", cv = FALSE,
                   empBayes = FALSE, preAssignment,
                   family = "binomial", prevfit,verbose=FALSE) {
  nvars <- ncol(mat)

  if(gamma < 0) gamma <- 0
  if(foreach::getDoParWorkers() == 1) {
    foreach::registerDoSEQ()
  }

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
    isingOffset[i] <- uniroot(f = function(off) mean(expit(eta + off)) - pResp,
                              interval = c(-50, 50))$root
  }
  # if(verbose) cat("\n")

  isingmat <- foreach(j = 1:ncol(mat), .combine = rbind) %dopar% {
    y <- as.vector(mat[, j])
    X <- as.matrix(mat[, -j])
    xcols <- colSums(X)
    if(family == "binomial") {
      regX <- X[, xcols >= 4]
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
                               intercept = FALSE)
      logliks <- 2 * (netfit$dev.ratio - 1) * netfit$nulldev
      dfs <- netfit$df
      ebic <- -logliks + dfs * log(nrow(mat) * (ncol(mat) - 1)^gamma)
      lambda <- netfit$lambda[which.min(ebic)]
    } else {
      netfit <- glmnet::cv.glmnet(regX, y, family = family, offset = off,
                                  intercept = FALSE)
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
      isingOffset[i] <- uniroot(f = function(off) mean(expit(eta + off)) - targetResp[i],
                                interval = c(-50, 50))$root
    }
    diag(isingmat) <- isingOffset
  }

  return(isingmat)
}
