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
  expect_is(flowReMixPFS(x = fit505, M = 5, stimVar = stimGroup, parentVar = parent),c("data.frame","tibble","data.table"))
  expect_true(all(c("parent","stimGroup","ptid","PFS")%in%colnames(flowReMixPFS(x = fit505, M = 5, stimVar = stimGroup, parentVar = parent))))
  expect_is(flowReMixPFS(x=fit505, M = 5, stimVar = stimGroup, parentVar = parent, parsefun = function(x,...)1),c("data.frame","tibble","data.table"))
  d1 = flowReMixPFS(x=fit505, M = 5, stimVar = stimGroup, parentVar = parent)
  d2 = flowReMixPFS(x=fit505, M = 5, stimVar = stimGroup, parentVar = parent, parsefun = function(x,...)5)
  expect_true(all(d1$PFS  <= d2$PFS))
  expect_true(all(table(weightForPFS(x=fit505,M = 5)) == c(85,27,8))) #distribution of weights
  expect_equal(degreeFromStringFun(c("env/23+/A+B+","env/A+B+","A+B+","A+","")),c("env/23+/A+B+"=2,"env/A+B+"=2,"A+B+"=2,"A+"=1,0))
  expect_error(flowReMixPFS(x=fit505, M = 5, stimVar = stimGroup, parentVar = par))
  expect_error(flowReMixPFS(x=fit505, M = 5, stimVar = stimG, parentVar = parent))
})
