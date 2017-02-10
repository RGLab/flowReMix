# Artificial data to work with ------------
n <- 250
nConditions <- 2
N <- sample(rv144$parentcount, n * nConditions) * 10
ptid <- as.vector(replicate(n, rep(runif(1), nConditions)))
condition <- as.vector(replicate(n, 1:nConditions - 1))
trueCondition <- condition
for(i in 1:length(trueCondition)) {
  if(condition[i] > 0) {
    trueCondition[i] <- 1 - rbinom(1, 1, 0.5)
  }
}
randomsig <- 0.000000000001
trueRandomsig <- randomsig
randomEffect <- as.vector(replicate(n, rep(rnorm(1, sd = randomsig), nConditions)))
M <- 10^3
trueM <- M
simdata <- data.frame(ptid, condition, N)
X <- model.matrix(~ factor(condition), simdata)
trueX <- model.matrix(~ factor(trueCondition), simdata)
beta <- c(-6, 0.5)
eta <- as.vector(trueX %*% beta) + randomEffect
expit <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log(x / (1 - x))
mu <- expit(eta)
p <- rbeta(length(mu), M * mu, M*(1 - mu))
boxplot(p ~ factor(expit(as.vector(X %*% beta))))
y <- rbinom(length(mu), N, p)
logpropmat <- matrix((y/N), ncol = nConditions, byrow = TRUE)
plot(logpropmat[, 1], logpropmat[, 2],
     ylim = c(min(logpropmat), max(logpropmat)),
     xlim = c(min(logpropmat), max(logpropmat)))
abline(a = 0, b = 1)
simdata$count <- y
waves <- as.vector(replicate(n, 1:nConditions))
simdata$waves <- waves
simdata$trueRandomEffect <- randomEffect
simdata$trueCondition <- trueCondition
simdata <- simdata[-1, ]
simdata$prop <- simdata$count / simdata$N
simdata$weights <- rep(1, nrow(simdata))
simdata$condition <- (simdata$condition)
simdata <- simdata[with(simdata, order(ptid, condition)), ]

compositeFit <- compositeMixture(count ~ condition,
                                     dformula = ~ 1,
                                     sub.population = NULL,
                                     N = N, id = ptid,
                                     waves = waves, data = simdata,
                                     treatment = condition,
                                     treatment.levels = 2,
                                     weights = NULL, init.phi = 1,
                                     init.beta = NULL, init.rho = 0, rho.fixed = FALSE,
                                     shrinkage.folds = 5,
                                     maxiter = 10, tol = 1e-02,
                                     computeT = TRUE,
                                     shrink.sandwich = TRUE)

c(M, (1 - compositeFit$rho[1])/ compositeFit$rho[1])
mean(trueCondition)*2
compositeFit$level.probabilities
compositeFit$coefTable

cummean <- function(x) cumsum(x) / 1:length(x)

trueCondition <- simdata$trueCondition
sprobMat <- cbind(compositeFit$subject.probs, trueCondition)[simdata$condition == 1, ]
sprobMat <- sprobMat[order(sprobMat[, 2]), ]
sprobMat <- cbind(sprobMat, 1 - cummean(sprobMat[,4]), cummean(sprobMat[, 2]))
round(sprobMat, 4)

sprobs <- fit$subject.probs
sprobs <- data.frame(round(data$count / data$parentcount, 4), sprobs, data$vaccine == "VACCINE")[seq(from = 2, to = nrow(sprobs), by = 2), ]
names(sprobs) <- c("prop", "subject", "nullprob", "altprob", "vaccine")
cor0 <- with(sprobs, cor(nullprob, vaccine))
cor1 <- with(sprobs, cor(altprob, vaccine))
if(cor0 > cor1) {
  sprobs <- sprobs[order(sprobs$nullprob), ]
  sprobs$nominalFDR <- round(cumsum(sprobs$nullprob) / 1:nrow(sprobs), 4)
  sprobs$respondProb <- sprobs$nullprob
} else {
  sprobs <- sprobs[order(sprobs$altprob), ]
  sprobs$nominalFDR <- round(cumsum(sprobs$altprob) / 1:nrow(sprobs), 4)
  sprobs$respondProb <- sprobs$altprob
}
sprobs$empFDR <- round(1 - cumsum(sprobs$vaccine) / 1:nrow(sprobs), 4)
sprobs$power <- with(sprobs, round(cumsum(vaccine)/sum(vaccine), 4))

require(pROC)
rocfit <- roc(sprobs$vaccine ~ sprobs$altprob)
par(mfrow = c(2, 1), mar = rep(2.5, 4))
plot(sprobs$nominalFDR, sprobs$empFDR,
     type= "l", xlim = c(0, 1), ylim = c(0, 1), col = "red")
lines(sprobs$nominalFDR, sprobs$power, col = "blue", lty = 2)
abline(a = 0, b = 1)
plot(rocfit, main = paste("AUC -", round(rocfit$auc, 4)))



# RV144 pvalue calibration ------------
### FULL RV144
data <- rv144
leaves <- unique(data$population)
data <- subset(data, stim != "sebctrl")
data$treatment <- as.numeric(data$stim == "env")
data <- subset(data, population %in% leaves[c(1, 2)])
data$ptid <- as.numeric(data$ptid)
data$ptid[data$vaccine == "VACCINE"] <- data$ptid[data$vaccine == "VACCINE"] * 10^4
data$prop <- data$count / data$parentcount
data$population <- as.factor(data$population)
data$stim[data$stim == "env"] <- "env"
data$stimB <- data$stim == "env"
data <- data[order(data$ptid, decreasing = FALSE), ]
fit <- compositeMixture(count ~ (treatment + age + gender),
                        N = parentcount,
                        data = data,
                        id = ptid,
                        shrink.sandwich = TRUE,
                        computeT = FALSE,
                        compute.offset = TRUE,
                        treatment = treatment,
                        treatment.levels = 2,
                        shrinkage.folds = 2,
                        sub.population = population,
                        maxiter = 100,
                        robust.param = 0.4,
                        tol = 10^-3)
fit$coefTable

geefit <- geepack::geeglm(count / parentcount ~ treatment + age + gender,
                                 data = data,
                                 id = ptid, family = "binomial")
summary(geefit)

truegeefit <- geepack::geeglm(count / parentcount ~ treatment*vaccine + age + gender,
                              data = data,
                              id = ptid, family = "binomial",
                              std.err = "j1s")

sprobs <- fit$subject.probs
sprobs <- data.frame(round(data$count / data$parentcount, 4), sprobs, data$vaccine == "VACCINE")[seq(from = 2, to = nrow(sprobs), by = 2), ]
names(sprobs) <- c("prop", "subject", "nullprob", "altprob", "vaccine")
sprobs$pvals <- round(fit$pvals[data$stim == "env"], 4)
sprobs$qvals <- round(p.adjust(sprobs$pvals, method = "fdr"), 4)
cor0 <- with(sprobs, cor(nullprob, vaccine, use = "complete.obs"))
cor1 <- with(sprobs, cor(altprob, vaccine, use = "complete.obs"))
if(fit$beta[2] > 0) {
  sprobsOrder <- order(sprobs$altprob, decreasing = TRUE)
  sprobs <- sprobs[sprobsOrder, ]
  sprobs$nominalFDR <- round(cumsum(sprobs$nullprob) / 1:nrow(sprobs), 4)
  sprobs$respondProb <- sprobs$nullprob
} else {
  sprobsOrder <- order(sprobs$nullprob, decreasing = TRUE)
  sprobs <- sprobs[sprobsOrder, ]
  sprobs$nominalFDR <- round(cumsum(sprobs$altprob) / 1:nrow(sprobs), 4)
  sprobs$respondProb <- sprobs$altprob
}
sprobs$empFDR <- round(1 - cumsum(sprobs$vaccine) / 1:nrow(sprobs), 4)
sprobs$power <- with(sprobs, round(cumsum(vaccine)/sum(vaccine), 4))

require(pROC)
par(mfrow = c(2, 2), mar = rep(2, 4))
colorIndex <- seq(from = 2, to = length(data$vaccine), by = 2)
props <- do.call("rbind", by(data$prop, data$ptid, function(x) x))
plot(log(props + 10^-5), col = (data$vaccine[colorIndex] == "VACCINE")+ 2)
abline(a = 0, b = 1)
rocfit <- roc(sprobs$vaccine ~ sprobs$altprob)
plot(sprobs$nominalFDR, sprobs$empFDR,
     type= "l", xlim = c(0, 1), ylim = c(0, 1), col = "red")
lines(sprobs$nominalFDR, sprobs$power, col = "blue", lty = 2)
abline(a = 0, b = 1)
plot(rocfit, main = paste("post AUC -", round(rocfit$auc, 4)))
rocfit <- roc(sprobs$vaccine ~ sprobs$pvals)
plot(rocfit, main = paste("pval AUC -", round(rocfit$auc, 4)))

props <- data.frame(props)
names(props) <- c("negctrl", "env")
props <- log(props + 10^-5)
props$vaccine <- as.vector(unlist(by(data$vaccine == "VACCINE", data$ptid, function(x) x[1])))
props <- props[sprobOrder, ]
props$posterior <- (sprobs$respondProb)
require(ggplot2)
ggplot(subset(props)) +
  geom_point(aes(x = negctrl, y = env, col = posterior, shape = !vaccine, size = !vaccine), fill = "white") +
  theme_bw() + geom_abline(intercept = 0, slope = 1) +
  scale_colour_gradientn(colours=rainbow(4))
  #facet_grid(vaccine ~ . , labeller = "label_both")

fit$coefTable
summary(geefit)
summary(truegeefit)



# Jackknifing
idlist <- unique(data$ptid)
betamat <- matrix(ncol = length(fit$beta), nrow = length(idlist))
for(i in 1:length(idlist)) {
  tempdat <- data[-which(data$ptid == idlist[[i]]), ]
  tempfit <- compositeMixture(count ~ treatment + stim + age + gender,
                          N = parentcount,
                          data = tempdat,
                          init.beta = fit$beta,
                          id = ptid,
                          shrink.sandwich = FALSE,
                          computeT = FALSE,
                          compute.offset = TRUE,
                          treatment = treatment,
                          treatment.levels = 2,
                          shrinkage.folds = 2,
                          sub.population = population,
                          maxiter = 50,
                          robust.param = 0.2,
                          tol = 10^-2)
  betamat[i, ] <- tempfit$beta
  print(c(i, i / length(idlist)))
  if(i > 1) {
    tempvar <- var(betamat[1:i, ]) * (nrow(betamat) - 1)
    print(tempfit$beta)
    print(sqrt(diag(tempvar)))
    print(fit$coefTable[, 2])
    print(fit$beta / sqrt(diag(tempvar)))
  }
}

