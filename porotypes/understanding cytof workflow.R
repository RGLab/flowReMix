generatePatient <- function(nConditions = 2) {
  subject_id <- runif(1)
  condition <- c(0:(nConditions - 1))
  sample_id <- interaction(subject_id, condition)
  y <- rbinom(nConditions, 1, 0.5)
  pdat <- data.frame(subject_id = subject_id, sample_id = sample_id, condition = condition, y = y)
  return(pdat)
}

library(lme4)
sampsize <- 50
data <- do.call("rbind", lapply(1:sampsize, function(x) generatePatient(2)))
fit <- glmer(y ~ condition + (1|subject_id) + (1|sample_id), data = data, family = binomial)

glmer(formula, data = NULL, family = gaussian, control = glmerControl(),
      start = NULL, verbose = 0L, nAGQ = 1L, subset, weights, na.action,
      offset, contrasts = NULL, mustart, etastart,
      devFunOnly = FALSE, …)
fit <- glmer
