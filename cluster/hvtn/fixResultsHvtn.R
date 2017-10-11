cpus <- 1
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)

library(flowReMix)

files <- as.list(dir(path = 'results', pattern="hvtn_*"))
dirfiles <- lapply(files, function(x) paste0('results/', x))
file <- dirfiles[[setting]]
load(file = file)

temp <- fit
temp$assignmentList <- NULL
temp$mhList <- NULL
temp$randomEffectSamp <- NULL
savefile <- paste("results2/", files[[setting]], sep ="")
save(temp, file = savefile)
file.remove(file)
rm(temp)

stability <- stabilityGraph(fit, type = "ising", cv = FALSE, reps = 200,
                            cpus = 1, AND = TRUE, gamma = 0.25)
stabfile <- paste("results2/stab_", files[[setting]], sep ="")
save(stability, file = stabfile)

