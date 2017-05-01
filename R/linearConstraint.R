linearMLEdensity <- function(constraint, y, product, contrastSig,
                             precision, productVar) {
  projMu <- y - (product - constraint) / (productVar) * contrastSig
  #mvtdens <- mvtnorm::dmvnorm(y, mean = projMu, sigma = sigma, log = TRUE)
  mvtdens <- mvtDens(y, projMu, precision)
  prob <- pnorm(threshold, mean = constraint, sd = sqrt(productVar), FALSE, TRUE)
  #print(c(constraint, mvtdens, prob))
  return(mvtdens - prob)
}

mvtDens <- function(y, mean, precision) {
  diff <- y - mean
  as.numeric(- 0.5 * t(diff) %*% precision %*% diff)
}

rejectSampleLinearConstraint <- function(mean, precision, contrast,
                                         threshold, df = 10^3,
                                         attempts = 10^4) {
  sample <- mvtnorm::rmvt(attempts, df = df, delta = mean, sigma = sigma)
  products <- sample %*% contrast
  keep <- products > threshold
  return(sample[keep, ])
}

linearConstraintMLE <- function(y, sigma, contrast, threshold, df = 10^3, alpha = 0.05,
                                FWcorrection = FALSE,
                                sampAttempts = 10^4 * 5) {
  p <- length(y)
  product <- as.numeric(t(y) %*% contrast)
  productVar <- t(contrast) %*% sigma %*% contrast
  contrastSig <- as.numeric(sigma %*% contrast)
  precision <- solve(sigma)

  mle <- optimize(f = linearMLEdensity, maximum = TRUE,
                  lower = 0, upper = product,
                  y = y, product = product, contrastSig = contrastSig,
                  precision = precision, productVar = productVar)$maximum

  mle <- y - (product - mle) / (productVar) * contrastSig
  mle <- pmax(mle, 0)

  sample <- rejectSampleLinearConstraint(rep(0, p), sigma, contrast, threshold,
                                         df = df,
                                         attempts = sampAttempts)
  constraints <- apply(sample, 1, function(z) {
    product <- t(z) %*% contrast
    optimize(f = linearMLEdensity, maximum = TRUE, lower = - 10 * product, upper = 10 * product,
             y = z, product = product, contrastSig = contrastSig,
             precision = precision, productVar = productVar)$maximum
    })
  sample <- t(sapply(constraints, function(mle) mle <- y - (product - mle) / (productVar) * contrastSig))
  sample <- apply(sample, 2, function(x) x - mean(x))

  globalPval <- pnorm(product, mean = 0, sd = sqrt(productVar), lower.tail = FALSE) /
    pnorm(threshold, mean = 0, sd = sqrt(productVar), lower.tail = FALSE)

  CI <- matrix(nrow = p, ncol = 2)
  naiveCI <- CI
  if(FWcorrection) {
    quantile <- qt(1 - alpha / p , df = df)
  } else {
    quantile <- qt(1 - alpha , df = df)
  }
  for(i in 1:p) {
    naiveCI[i, ] <- y[i] + c(-1, 1) * sqrt(sigma[i, i]) * quantile
  }

  if(globalPval < alpha^2) {
    CI <- naiveCI
  } else {
    if(FWcorrection) {
      qalpha <- (alpha - alpha^2) / p
    } else {
      qalpha <- (alpha - alpha^2)
    }
    for(i in 1:p) {
      CI[i, ] <- mle[i] - quantile(sample[, i], c(1 - qalpha / 2, qalpha / 2))
    }
  }


  return(list(naive = y, naiveCI = naiveCI, mle = mle, CI = CI))
}

# p <- 5
# rho <- 0.5
# sigma <- matrix(rho, ncol = p, nrow = p)
# diag(sigma) <- 1
# mu <- abs(rnorm(p, sd = 0.5))
# contrast <- rep(1, p) / p
# var <- as.numeric(t(contrast) %*% sigma %*% contrast)
# threshold <- 1.64 * sqrt(var)
# y <- as.numeric(mvtnorm::rmvnorm(1, mean = mu, sigma = sigma))
# obs <- as.numeric(contrast %*% y)
# while(obs < threshold) {
#   y <- as.numeric(mvtnorm::rmvnorm(1, mean = mu, sigma = sigma))
#   obs <- contrast %*% y
# }
