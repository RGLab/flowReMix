library(flowReMix)
ncpus <- 15

filenames <- as.list(dir(path = 'results', pattern="rv144_1_*"))
filenames <- lapply(filenames, function(x) paste0('results/', x))

post <- list()
random <- list()
assign <-list()
for(i in 1:length(filenames)) {
  load(file = filenames[[i]])
  post[[i]] <- fit$posteriors[, -1]
  random[[i]] <- fit$randomEffectSamp
  assign[[i]] <- fit$assignmentList
}
post <- Reduce("+", post) / length(filenames)
assign <- do.call("c", assign)
random <- do.call("c", random)
fit$posteriors[, -1] <- post
fit$randomEffectSamp <- random
fit$assignmentList <- assign

stability <- stabilityGraph(fit, type = "ising", cv = TRUE, reps = 1000,
                            cpus = ncpus, AND = TRUE)
save(stability, file = "results/rv144AggreageStability1.Robj")

stability <- stabilityGraph(fit, type = "randomEffects", cv = TRUE, reps = 1000,
                            cpus = ncpus, AND = TRUE)
save(stability, file = "results/rv144AggreageRandom1.Robj")
