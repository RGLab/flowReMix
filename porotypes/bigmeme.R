nsubjects <- 236
npost <- 3
nsamp <- 100
keepEach <- 5
nSubsets <- 240
niter <- 50
nrows <- nsubjects * npost * nsamp * niter / keepEach
bigmat <- big.matrix(init = 0, nrow = nrows, ncol = nSubsets)
GetMatrixSize(bigmat) / 10^6
