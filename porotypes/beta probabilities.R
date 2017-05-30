range <- seq(from = 0, to = 1, length.out = 10000)
alpha <- 1
beta <- 10
par(mfrow = c(1, 1))
plot(range, dbeta(range, alpha, beta), type = "l")

msize <- 0:100
modelprobs <- beta(alpha + msize, beta + max(msize) - msize) / beta(alpha, beta)
modelprobs <- 1/(1 + 3*(0:100))
modelprobs <- modelprobs / sum(modelprobs)
plot(diff(log(modelprobs)), type = "l")

