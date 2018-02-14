context("scores")
library(flowReMix)
library(testthat)
data(fit505)
test_that("flowReMix objects are detected correctly",{
  expect_error(flowReMix:::.isFlowRemix(NULL))
  expect_equal(NULL,flowReMix:::.isFlowRemix(fit505))
})

test_that("getIsing, getPosteriors, getSubsets work correctly",{
  expect_equal(getSubsets(x = fit505),colnames(fit505$posteriors)[-1L])
  expect_equal(getPosteriors(x=fit505),fit505$posteriors)
  expect_equal(getIsing(x=fit505),fit505$isingStability)
})


test_that("flowReMix object has required fields",{
  expect_true(all(c("coefficients","posteriors","data","call","isingStability","levelProbs","subject_id","control")%in%names(fit505)))
})

test_that("Polyfunctionality Score Functions work correctly", {
  expect_is(flowReMixPFS(x = fit505, M = 5, stimVar = stimGroup, subsetVar=subset, split="\\+",parentVar = parent),c("data.frame","tibble","data.table"))
  expect_true(all(c("parent","stimGroup","ptid","PFS")%in%colnames(flowReMixPFS(x = fit505,  subsetVar=subset,split="\\+",M = 5, stimVar = stimGroup, parentVar = parent))))
  expect_warning(flowReMixPFS(x=fit505, M = 5,  split="\\+",stimVar = stimGroup, subsetVar=subset, parentVar = parent, parsefun = function(x,...)1))
  d1 = flowReMixPFS(x=fit505, M = 5,  split="\\+",stimVar = stimGroup, subsetVar=subset, parentVar = parent)
  suppressWarnings({d2 = flowReMixPFS(x=fit505, M = 5,  split="\\+",subsetVar=subset,stimVar = stimGroup, parentVar = parent, parsefun = function(x,...)5)})
  expect_true(all(d1$PFS  <= d2$PFS))
  expect_true(all(table(weightForPFS(x=fit505,M = 5,split = "\\+")) == c(85,27,8))) #distribution of weights
  expect_equal(degreeFromStringFun(c("env/23+/A+B+","env/A+B+","A+B+","A+","")),c("env/23+/A+B+"=2,"env/A+B+"=2,"A+B+"=2,"A+"=1,0))
  expect_error(flowReMixPFS(x=fit505, subsetVar=subset,M = 5,split="\\+", stimVar = stimGroup, parentVar = par))
  expect_error(flowReMixPFS(x=fit505, subsetVar=subset, M = 5, split="\\+",stimVar = stimG, parentVar = parent))
})

test_that("scatterplots work and accept subset argument",{
  data(fit505)
  subs = getSubsets(fit505)[1:3]
  expect_is({testplot = plot(fit505,subsets=subs,target=vaccine,type="scatter")},"ggplot")
  expect_equal(nlevels(factor(testplot$data$sub.population)),3)
})

test_that("RV144 polyfunctionality scores are reproducible",{
  data("flowremix_testdata")
  data("rv144_pfs_test")
  data("fitcontrol")
  control$updateLag=4
  control$ncores=4
  control$threads=128
  fit <- flowReMix(cbind(count, parentcount - count) ~ treatment,
                   subject_id = ptid,
                   cell_type = subset,
                   cluster_variable = treatment,
                   data = booldata,
                   covariance = "sparse",
                   ising_model = "sparse",
                   regression_method = "robust",
                   iterations = 10,
                   cluster_assignment = TRUE,
                   parallel = TRUE, keepSamples = TRUE,
                   verbose = FALSE, control = control,
                   newSampler = FALSE);
  fit$data$parent="CD4";
  suppressWarnings({results = flowReMixPFS(fit,parentVar = "parent",subsetVar=subset, stimVar = "stim", split=",",M=6,outcomeVar = "hiv")%>%dplyr::inner_join(pfs_test%>%dplyr::select(PFS_COMPASS=PFS,ptid=PTID)%>%dplyr::mutate(ptid=factor(ptid)))})
  expect_gt(cor(results$PFS,results$PFS_COMPASS),0.97)
})
