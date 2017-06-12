BBreg <- function(data, formula, weights) {
  # Prototype data -----
  # n <- 100
  # p <- 10
  # X <- matrix(rnorm(n * p), nrow = n)
  # X <- cbind(1, X)
  # beta <- c(-3, runif(p, min = -0.2, max = 0.2))
  # offset <- rnorm(n, sd = 0.2)
  # eta <- as.vector(X %*% beta) + offset
  # mu <- 1 / (1 + exp(-eta))
  # M <- 10^2
  # prop <- rbeta(n, M * mu, M * (1 - mu))
  # N <- rpois(n, rnorm(n, 10^4, 10^3))
  # y <- rbinom(n, N, prop)
  # data <- data.frame(y, N, offset, V = X)
  # formula <- as.formula(paste("cbind(y, N - y) ~ offset(offset)",
  #                             paste("+ V.", 1:ncol(X), collapse = "", sep = "")))
  # weights <- runif(n)

  # Fitting the model
  glmfit <- glm(formula, data = data, family = "binomial", weights = weights)
  coef <- coef(glmfit)
  toremove <- is.na(coef)
  coef <- coef[!toremove]
  X <- model.matrix(formula, data)
  X <- X[, !toremove]
  response <- model.response(glmfit$model)
  weights <- model.weights(glmfit$model)
  offset <- model.offset(glmfit$model)
  y <- response[, 1]
  N <- y + response[, 2]
  mu <- predict(glmfit, type = "response")

  M <- dispersionMLE(y, N, mu)
  param <- c(coef, M)
  param <- optim(param, fn = evalBB, gr = compgrad,
        method = "L-BFGS-B",
        lower = c(rep(-10^5, ncol(X)), 10),
        y = y, N = N, X = X, offset = offset,
        weights = weights)
  coef <- param$par[-length(param$par)]
  eta <- X %*% coef
  coef <- coef(glmfit)
  coef[!is.na(coef)] <- param$par[1:ncol(X)]
  coef[is.na(coef)] <- 0
  M <- param$par[length(param$par)]

  result <- glmfit
  result$M <- M
  result$coefficients <- coef
  result$linear.predictors <- eta
  class(result) <- c("bbreg", class(result))
  return(result)
}

compgrad <- function(param, y, N, X, offset, weights) {
  M <- param[length(param)]
  coef <- param[1:(length(param) -1)]
  # coef <- param
  eta <- as.numeric(X %*% coef) + offset
  mu <- 1 / (1 + exp(- eta))
  diff <- N - y
  alpha <- M * mu
  beta <- M * (1 - mu)
  digam <- digamma(y + alpha) - digamma(diff + beta) - digamma(alpha) + digamma(beta)
  digam <- digam * M
  Mgrad <- digamma(y + alpha) * mu + digamma(N - y + beta) * (1 - mu) - digamma(N + M)
  Mgrad <- Mgrad + digamma(M) - digamma(alpha) * mu - digamma(beta) * (1 - mu)
  diexp <- exp(eta) / (exp(eta) + 1) ^ 2
  grad <- t(X) %*% (diexp * digam * weights)
  Mgrad <- sum(Mgrad * weights)

  return(-c(grad, Mgrad))
  # return(- grad)
}

evalBB <- function(param, y, N, X, offset, weights) {
  M <- param[length(param)]
  coef <- param[1:(length(param) -1)]
  # coef <- param
  eta <- as.numeric(X %*% coef) + offset
  mu <- 1 / (1 + exp(- eta))
  return(-sum(vecBetaBinomDens(y, N, mu, M) * weights))
}

predict.bbreg <- function(object) {

}

