# Adapted from: http://www.math.yorku.ca/~hkj/Software/gren.R
findGren  <-  function(p) {
  n <- length(p)
  pts <- cbind(c(0:n, n), c(0, cumsum(p), 0))

  hpts <- chull(pts)
  hpts <- c(hpts, hpts[1])

  hpairs <- matrix(0, length(hpts) - 1, 2)
  for(i in 1:(length(hpts)-1)){
    hpairs[i, ]	<-	c(hpts[i + 1], hpts[i])
  }
  m <- length(pts[, 1])
  hpairs <- hpairs[which(hpairs[, 1] != m), ]
  hpairs  <-  hpairs[which(hpairs[, 2] != m), ]
  if(length(hpairs) == 2){
    hpairs <-	matrix(hpairs, 1, 2)
  } else 	{
    hpairs <-	hpairs[order(hpairs[, 1]), ]
  }

  s.num <- pts[hpairs[, 1], 2] - pts[hpairs[, 2], 2]
  s.denom <- pts[hpairs[, 1], 1] - pts[hpairs[, 2], 1]
  slopes <- s.num / s.denom
  h.hat <- numeric()
  for(i in 1:length(slopes)){
    h.hat	<- c(h.hat, rep(slopes[i], hpairs[i, 1] - hpairs[i, 2]))
  }

  return(h.hat)
}

estimateMonotoneProbs <- function(samp, method = "arrange") {
  nSubsets <- ncol(samp)
  samp <- rowSums(samp)
  countTable <- sapply(0:nSubsets, function(x) mean(x == samp))
  if(method == "gren") {
    modelprobs <- findGren(countTable)
  } else if(method == "arrange") {
    modelprobs <- sort(countTable, decreasing = TRUE)
  }

  modelprobs <- pmax(modelprobs, min(modelprobs[modelprobs > 0]) / 2)
  modelprobs <- modelprobs / sum(modelprobs)
  return(modelprobs)
}
