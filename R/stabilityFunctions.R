#' @export
stabilityGraph <- function(obj, type = c("ising", "randomEffects"),
                           cv = FALSE, reps = 100, cpus = 1,
                           gamma = 0.9, AND = TRUE) {
  type <- type[1]

  if(type == "ising") {
    samples <- obj$assignmentList
    family <- "binomial"
  } else if(type == "randomEffects") {
    samples <- obj$randomEffectSamp
    family <- "gaussian"
  } else {
    stop("Method not supported!")
  }
  names(samples) <- sapply(names(samples), function(x) strsplit(x, "%%%")[[1]][[1]])
  samples <- lapply(unique(names(samples)), function(x) {
    do.call("rbind", samples[names(samples) == x])
  })
  subsets <- names(obj$coefficients)

  nsubsets <- ncol(samples[[1]])
  countCovar <- matrix(0, nrow = nsubsets, ncol = nsubsets)
  if(cpus == 1) {
    foreach::registerDoSEQ()
  } else {
    doParallel::registerDoParallel(cores = cpus)
  }

  perc <- 0.1
  cat("Progress: ")
  for(i in 1:reps) {
    mat <- t(sapply(samples, function(x) x[sample(1:nrow(x), 1), ]))
    colnames(mat) <- subsets
    coefs <- raIsing(mat, AND = AND, gamma = gamma, family = family,
                     method = "sparse", cv = cv)
    countCovar <- countCovar + (coefs != 0) * sign(coefs)
    if(i / reps > perc & perc < 1) {
      cat(perc * 100, "% ", sep = "")
      perc <- perc + 0.1
    }
  }
  cat("100% \n")
  doParallel::stopImplicitCluster()

  props <- countCovar / reps
  colnames(props) <- names(obj$coefficients)
  rownames(props) <- colnames(props)
  diag(props) <- 0
  results <- list()
  results$network <- props
  results$frequencies <- table(props)
  results$type <- type
  class(results) <- c("flowReMix_stability")
  return(results)
}

#' @export
plot.flowReMix_stability <- function(obj, threshold = 0.5, plotAll = FALSE,
                                     fill = NULL, fillRange = NULL, fillPalette = NULL,
                                     title = TRUE) {
  measure <- fill
  props <- obj$network
  props[abs(props) < threshold] <- 0
  network <- props
  if(!plotAll) {
    keep <- apply(network, 1, function(x) any(abs(x) >= threshold))
  }
  network <- network[keep, keep]
  net <- network::network(props)
  subsets <- colnames(props)
  nodes <- GGally::ggnet2(network, label = subsets[keep])$data
  edges <- matrix(nrow = sum(network != 0)/2, ncol = 5)

  p <- nrow(network)
  row <- 1
  for(j in 2:p) {
    for(i in 1:(j-1)) {
      if(network[i, j] != 0) {
        edges[row, ] <- unlist(c(nodes[i, 6:7], nodes[j, 6:7], network[i, j]))
        row <- row + 1
      }
    }
  }

  edges <- data.frame(edges)
  names(edges) <- c("xstart", "ystart", "xend", "yend", "width")
  if(!is.null(measure)) {
    nodes$measure <- measure[keep]
  }

  names(edges)[5] <- "Dependence"
  edges$Correlation <- c("positive")
  edges$Correlation[edges$Dependence < 0] <- "negative"
  edges$Correlation <- factor(edges$Correlation, levels = c("positive", "negative"))
  edgecol <- scale_colour_manual(name = "Correlation", values = c("dark green", "dark red"))
  edges$Dependence <- abs(edges$Dependence)

  library(ggplot2)
  figure <- ggplot() +
    edgecol +
    geom_segment(data = edges, aes(x = xstart, y = ystart,
                                   xend = xend, yend = yend,
                                   col = Correlation,
                                   alpha = Dependence),
                 size = 1) +
    scale_alpha_continuous(name = "Edge Propensity")
    if(!is.null(measure)) {
      measureName <- as.character(match.call()$fill)
      measureName <- measureName[length(measureName)]
      if(is.null(fillRange)) {
        fillRange <- c(min(measure), max(measure))
      }
      if(is.null(fillPalette)) {
        fillPalette <- rainbow(4)
      }

      figure <- figure + geom_point(data = nodes, aes(x = x, y = y, fill = nodes$measure), shape = 21, size = 8, col = "grey") +
        scale_fill_gradientn(limits = fillRange, colours = fillPalette, name = measureName)
    } else {
      figure <- figure + geom_point(data = nodes, aes(x = x, y = y), shape = 21, size = 8, col = "grey",
                                    fill = "sky blue")
    }
  figure <- figure + scale_shape(solid = FALSE) +
    geom_text(data = nodes, aes(x = x, y = y, label = nodes$label), size = 1.8) +
    theme(axis.line=element_blank(), axis.text.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(), plot.background=element_blank())

  if(is.character(title)) {
    figure <- figure + ggtitle(title)
  } else if(obj$type == "ising" & title) {
    figure <- figure + ggtitle("Estimated Dependence Structure for Ising Model")
  } else if(obj$type == "randomEffects" & title) {
    figure <- figure + ggtitle("Estimated Dependence Structure for Random Effects")
  }

  return(figure)
}
