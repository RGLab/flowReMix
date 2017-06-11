BBreg <- function(data, formula) {
  # Prototype data -----
  # n <- 100
  # p <- 10
  # X <- matrix(rnorm(n * p), nrow = n)
  # X <- cbind(1, X)
  # beta <- c(-3, runif(p, min = -0.2, max = 0.2))
  # offset <- rnorm(n, sd = 0.2)
  # eta <- as.vector(X %*% beta) + offset
  # mu <- 1 / (1 + exp(-eta))
  # M <- 10^3
  # prop <- rbeta(n, M * mu, M * (1 - mu))
  # N <- rpois(n, rnorm(n, 10^4, 10^3))
  # y <- rbinom(n, N, prop)
  # data <- data.frame(y, N, offset, V = X)
  # formula <- as.formula(paste("cbind(y, N - y) ~ offset(offset)",
  #                             paste("+ V.", 1:ncol(X), collapse = "", sep = "")))

  # Fitting the model
  glmfit <- glm(formula, data = data, family = "binomial")
  coef <- coef(glmfit)
  toremove <- is.na(coef)
  coef <- coef[!toremove]
  X <- model.matrix(formula, data)
  X <- X[, !toremove]
  response <- model.response(glmfit$model)
  y <- response[, 1]
  N <- y + response[, 2]
  mu <- predict(glmfit, type = "response")

  M <- dispersionMLE(y, N, mu)
  param <- c(coef, M)
  param <- optim(param, fn = evalBB, gr = compgrad,
        method = "L-BFGS-B",
        lower = c(rep(-10^5, ncol(X)), 50),
        y = y, N = N, X = X, offset = offset)
  mu <- as.numeric(1 / (1 + exp(-(X %*% param$par))))
  M <- dispersionMLE(y, N, mu)
  cbind(coef, beta, param$par[-length(param$par)])
  coef <- coef(glmfit)
  coef[!is.na(coef)] <- param$par[1:ncol(X)]
  coef[is.na(coef)] <- 0
  M <- param$par[length(param$par)]
  return(list(coef = coef, M = M))
}

compgrad <- function(param, y, N, X, offset) {
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
  grad <- t(X) %*% (diexp * digam)
  Mgrad <- sum(Mgrad)

  return(-c(grad, Mgrad))
  # return(- grad)
}

evalBB <- function(param, y, N, X, offset) {
  M <- param[length(param)]
  coef <- param[1:(length(param) -1)]
  # coef <- param
  eta <- as.numeric(X %*% coef) + offset
  mu <- 1 / (1 + exp(- eta))
  -sum(vecBetaBinomDens(y, N, mu, M))
}

