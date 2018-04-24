library(flowReMix)
library(magrittr)
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)

# loading file
filenames <- as.list(dir(path = 'results', pattern="HVTN105_C_110subsetsTEST_iterations30disp*"))
filenames <- lapply(filenames, function(x) paste0('results/', x))
filenames <- as.character(filenames)
fit <- readRDS(filenames[setting])
fit$isingStability <- stabilityGraph(fit, type = "ising", cv = TRUE,
                                     reps = 100, cpus = 1)
fit$assignmentList <- NULL
fit$randomEffectSamp <- NULL
fit$data <- NULL
file <- filenames[setting]
file <- strsplit(file, "/", fixed = TRUE)[[1]]
file <- file[length(file)]
file <- paste("results2", file, sep = "/")
saveRDS(fit, file = file)


