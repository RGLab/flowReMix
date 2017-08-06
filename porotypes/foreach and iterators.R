system.time(l <- lapply(1:10^4, function(x) rnorm(100)))
liter <- iterators::iter(l)
foreachResult <- foreach::foreach(obj = liter, .combine = "c") %dopar% {
  a <- mean(obj)
  return(a)
}

cbind(foreachResult, sapply(l, mean))

system.time(l <- lapply(1:10^5, function(x) rnorm(100)))
system.time(l2 <- lapply(1:10^5, function(x) rnorm(100)))
system.time(l3 <- mapply(function(x,y) list(x, y), l, l2, SIMPLIFY = FALSE))
