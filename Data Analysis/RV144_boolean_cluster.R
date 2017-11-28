library(flowReMix)
library(foreach)
library(doParallel)
library(doRNG)

assign <- function(x) {
  x$prop <- x$count / x$parentcount
  assign <- as.numeric(by(x, x$subset, function(y) y$prop[1] > y$prop[2]))
  assign[assign == 1] <- -1
  result <- data.frame(ptid = x$ptid[1], subset = unique(x$subset), assign = assign)
  return(result)
}

require(pROC)
require(reshape2)
load("../data/rv144_booleans.rda")
bySubset <- by(data.frame(booleans$stim, booleans$nonstim), booleans$Subset, function(x) x)
largerThanThershold <- sapply(bySubset, function(x) colSums(x >5))

booldata <- melt(booleans, c("PTID", "Subset"))
names(booldata)[3:4] <- c("stim", "count")

booldata <- by(booldata, INDICES = list(booldata$PTID, booldata$stim), function(x) {
  x$parentcount <- sum(x$count)
  return(x)
})
booldata <- do.call("rbind", booldata)

booldata <- subset(booldata, Subset != "!TNFa&!IFNg&!IL4&!IL2&!CD154&!IL17a")
booldata$treatment <- as.numeric(booldata$stim == "stim")
uniquepop <- unique(booldata$Subset)
booldata <- with(booldata, booldata[order(Subset, PTID, stim, decreasing = FALSE), ])
booldata <- subset(booldata, !is.na(Subset))
allsubset <- booldata
booldata <- with(booldata, booldata[order(Subset, PTID, stim, decreasing = FALSE), ])

# Naming ------------------
subsets <- unique(booldata$Subset)
booldata$Subset <- as.character(booldata$Subset)
nfunctions <- numeric(length(subsets))
for(i in 1:length(subsets)) {
  split <- strsplit(as.character(subsets[i]), "&")[[1]]
  first <- substr(split, 1, 1)
  nfunction <- sum(first != "!")
  nfunctions[i] <- nfunction
  name <- paste(split[first != "!"], collapse = ",")
  booldata$nfunction[booldata$Subset == subsets[[i]]] <- nfunction
  booldata$Subset[booldata$Subset == subsets[[i]]] <- name
}
subsets <- unique(booldata$Subset)
booldata <- with(booldata, booldata[order(Subset, PTID, stim, decreasing = FALSE), ])
names(booldata) <- tolower(names(booldata))


# Getting vaccine information --------------------
data("rv144")
rv144 <- rv144[order(rv144$ptid), ]
vaccine <- (by(rv144, rv144$ptid, function(x) x$vaccine[1] == "VACCINE"))
vaccine <- data.frame(ptid = names(vaccine), vaccine = as.numeric(vaccine))
vaccinemat <- vaccine[vaccine$ptid %in% booldata$ptid, ]

# Getting infection status
data("rv144_correlates_data")
correlates <- rv144_correlates_data
correlates <- correlates[order(as.character(correlates$PTID)), ]
infection <- correlates$infect.y
names(correlates)[1] <- "ptid"
vaccine <- correlates[, c(1, 62, 64)]
vaccine$ptid <- factor(vaccine$ptid, levels = levels(booldata$ptid))
booldata <- with(booldata, booldata[order(subset, ptid, stim, decreasing = FALSE), ])
booldata <- merge(booldata, vaccine, sort = FALSE)
booldata$hiv <- NA
booldata$hiv[booldata$infect.y == "INFECTED"] <- TRUE
booldata$hiv[booldata$infect.y == "NON-INFECTED"] <- FALSE
booldata <- with(booldata, booldata[order(subset, ptid, stim, decreasing = FALSE), ])

# Analysis -------------
library(flowReMix)

pargrid = expand.grid(niter=c(20,50,100,500,1000),nsamp = c(20,50,100),prior=c(0,1,2,4))
setting = Sys.getenv('SLURM_ARRAY_TASK_ID')
cores = Sys.getenv('SLURM_CPUS_PER_TASK')
config = pargrid[setting,]
npost <- 1
niter  = as.integer(config[1])
prior = as.numeric(config[3])
nsamp = as.integer(config[2])
 control <- flowReMix_control(updateLag = 10, nsamp = nsamp, initMHcoef = 2.5,
                              keepEach = 1,
                              nPosteriors = npost, centerCovariance = TRUE,
                              maxDispersion = 1000, minDispersion = 10^7,
                              randomAssignProb = 10^-8, intSampSize = 50,
                              lastSample = 4, isingInit = -log(99),
                              ncores = cores,
                              preAssignCoefs = 0,
                              prior = prior, isingWprior = TRUE,
                              markovChainEM = TRUE,
                              initMethod = "robust",
                              zeroPosteriorProbs = TRUE)

 booldata$subset <- factor(booldata$subset)
 preAssignment <- do.call("rbind", by(booldata, booldata$ptid, assign))
 booldata$stim <- factor(booldata$stim, levels = c("nonstim", "stim"))
 system.time(fit <- flowReMix(cbind(count, parentcount - count) ~ stim,
                              subject_id = ptid,
                              cell_type = subset,
                              cluster_variable = stim,
                              data = booldata,
                              covariance = "sparse",
                              ising_model = "sparse",
                              regression_method = "robust",
                              iterations =  niter,
                              cluster_assignment = TRUE,
                              parallel = TRUE,
                              verbose = FALSE, control = control))
 saveRDS(fit,paste0("output/rv144_cluster_niter_",niter,"_nsamp_",nsamp,"_prior_",prior,".rds"))
# plot(fit, type = "scatter", target = vaccine)
# plot(fit, type = "ROC", target = vaccine, direction = "<",varname="vaccine")
# plot(fit, type = "ROC", target = hiv, direction = ">",varname="hiv")
# summary(fit,target=hiv)%>%arrange(-auc)%>%head(10)
# summary(fit,target=vaccine)%>%arrange(-auc)%>%head(10)
# M <- 6
# weights <- nfunctions / (choose(M, nfunctions))
# poly <- apply(fit$posteriors[, -1], 1, function(x) weighted.mean(x, weights))
# post = gather(data=fit$posteriors,key = subset,value = posterior,-1)%>%left_join(preAssignment)
#
# post%>%mutate(posterior_new = ifelse(assign==0,0,posterior))%>%group_by(ptid,subset)%>%mutate(weight=length(stringr::str_split(pattern = ",",subset)[[1]]))%>%mutate(weight = weight/choose(6,weight))%>%group_by(ptid)%>%
#   summarize(poly_postassign=weighted.mean(posterior_new,weights),poly=weighted.mean(posterior,weights))%>%left_join(correlates%>%select(ptid,PFS))%>%
#   do(data.frame(correlation=cor(.$poly,.$PFS)))
# post%>%mutate(posterior_new = ifelse(assign==0,0,posterior))%>%group_by(ptid,subset)%>%mutate(weight=length(stringr::str_split(pattern = ",",subset)[[1]]))%>%mutate(weight = weight/choose(6,weight))%>%group_by(ptid)%>%
#   summarize(poly_postassign=weighted.mean(posterior_new,weights),poly=weighted.mean(posterior,weights))%>%left_join(correlates%>%select(ptid,PFS))%>%
#   ggplot()+geom_point()+aes(x=PFS,y=poly)+geom_smooth(method="lm")+theme_classic()+scale_y_continuous("Total Response")
#
# #outliers
# post%>%mutate(posterior_new = ifelse(assign==0,0,posterior))%>%group_by(ptid,subset)%>%mutate(weight=length(stringr::str_split(pattern = ",",subset)[[1]]))%>%mutate(weight = weight/choose(6,weight))%>%group_by(ptid)%>%
#   summarize(poly_postassign=weighted.mean(posterior_new,weights),poly=weighted.mean(posterior,weights))%>%left_join(correlates%>%select(ptid,PFS))%>%mutate(d=as.vector(scale(PFS)-scale(poly_postassign))^2)%>%ungroup%>%dplyr::arrange(-d)%>%select(ptid)%>%head(5)
# fit$posteriors["P1007",]%>%gather(subset,posterior,-1)%>%
# left_join(fit$data%>%filter(ptid=="P1007"))%>%select(ptid,subset,posterior,stim,count,parentcount)%>%mutate(prop=count/parentcount)%>%select(ptid,subset,posterior,prop,stim)%>%spread(stim,prop)%>%mutate(diff=stim-nonstim)%>%filter(diff<=0)
#
# wts = unlist(lapply(strsplit(colnames(fit$posteriors)[-1],","),length))/choose(6,unlist(lapply(strsplit(colnames(fit$posteriors)[-1],","),length)))
# plot(fit,type="boxplot",group="all",target=hiv,jitter=TRUE,test="wilcoxon",one_sided=TRUE,weights=list(wts))
#
#
# summary(glm(infect.x ~ scale(poly) + IgAprim + risk.medium + risk.high + sex,
#             family = "binomial",
#             data = post%>%mutate(posterior_new = ifelse(assign==0,0,posterior))%>%group_by(ptid,subset)%>%mutate(weight=length(stringr::str_split(pattern = ",",subset)[[1]]))%>%mutate(weight = weight/choose(6,weight))%>%group_by(ptid)%>%
#               summarize(poly_postassign=weighted.mean(posterior_new,weights),poly=weighted.mean(posterior,weights))%>%left_join(correlates%>%select(ptid,PFS,infect.x,IgAprim,sex,risk.high,risk.medium))))
#
