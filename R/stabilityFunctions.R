#' @name stabilityGraph
#' @title compute a stability graph for the Ising model from a flowReMix fit.
#' @import ggplot2
#' @importFrom network network
#' @importFrom GGally ggnet2
#' @importFrom viridis plasma
#' @importFrom igraph graph.adjacency
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
    registerDoRNG()
  }
  set.seed(100)

  # perc <- 0.1
  # requireNamespace("progress")
  # pb = progress_bar$new(total=reps);
  cluster_res = foreach(i = 1:reps) %dorng% {
    mat <- t(sapply(samples, function(x) x[sample(1:nrow(x), 1), ,drop=FALSE]))
    colnames(mat) <- subsets
    coefs <- raIsing(mat, AND = AND, gamma = gamma, family = family,
                     method = "sparse", cv = cv)
    countCovar <- countCovar + (coefs != 0) * sign(coefs)
    # if(i / reps > perc & perc < 1) {
    #   cat(perc * 100, "% ", sep = "")
    #   perc <- perc + 0.1
    # }
    return(countCovar)
  }
  countCovar = Reduce(x = cluster_res, f=function(x,y)x+y)
  # cat("100% \n")
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
plotRawGraph <- function(obj, graph = c("ising"), threshold = 0.5, plotAll = FALSE,
                         fill = NULL, fillName = NULL,
                         fillRange = NULL, fillPalette = NULL,
                         title = TRUE, normalize = FALSE,
                         count = TRUE, label_size = 1.8,seed=100) {
  if(graph == "ising") {
    ising <- obj$isingAvg
    if(count) {
      isingZ <- obj$isingCount
    } else if(normalize) {
      isingZ <- ising / sqrt(obj$isingVar)
    } else {
      isingZ <- ising
    }
  } else if(graph == "randomEffects") {
    ising <- obj$invCovAvg
    isingZ <- ising / sqrt(obj$invCovVar)
    if(normalize) {
      isingZ <- ising / sqrt(obj$isingVar)
    } else {
      isingZ <- ising
    }
  } else {
    stop("Graph type not supported!")
  }
  isingZ[is.nan(isingZ)] <- 0
  diag(isingZ) <- 0
  edges <- isingZ[lower.tri(isingZ)]
  nonZero <- abs(edges) > 0
  edges[nonZero] <- sapply(edges[nonZero], function(x) mean(abs(edges[nonZero]) <= abs(x)))
  isingZ[lower.tri(isingZ)] <- edges
  for(i in 1:(nrow(isingZ) - 1)) {
    for(j in (i+1):ncol(isingZ)) {
      isingZ[i, j] <- isingZ[j, i]
    }
  }
  network <- list()
  rownames(isingZ) <- names(fit$coefficients)
  colnames(isingZ) <- names(fit$coefficients)
  network$network <- isingZ
  network$type <- "ising"

  if(is.null(fillName) & !is.null(fill)) {
    fillName <- as.character(match.call()$fill)
    fillName <- fillName[length(fillName)]
  }

  network$raw <- TRUE

  figure <- plot.flowReMix_stability(network, threshold = threshold, plotAll = plotAll,
                                     fill = fill, fillRange = fillRange, fillPalette = fillPalette,
                                     fillName = fillName,
                                     title = title, label_size = label_size, seed=seed)
  return(figure)
}

#' @export
plot.flowReMix_stability <- function(obj, threshold = 0.5, plotAll = FALSE,
                                     fill = NULL, fillName = NULL,
                                     fillRange = NULL, fillPalette = NULL,
                                     title = TRUE, label_size = 1.8,seed=100) {
  set.seed(seed)
  require(ggplot2)
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

  if(is.null(obj$raw)) {
    alphaTitle <- "Edge Propensity"
  } else {
    alphaTitle <- "Edge Quantile"
  }

  library(ggplot2)
  figure <- ggplot() +
    edgecol +
    geom_segment(data = edges, aes(x = xstart, y = ystart,
                                   xend = xend, yend = yend,
                                   col = Correlation,
                                   alpha = Dependence),
                 size = 1) +
    scale_alpha_continuous(name = alphaTitle)
    if(!is.null(measure)) {
      if(is.null(fillName)) {
        measureName <- as.character(match.call()$fill)
        measureName <- measureName[length(measureName)]
      } else {
        measureName <- fillName
      }
      if(is.null(fillRange)) {
        fillRange <- c(min(measure), max(measure))
      }
      if(is.null(fillPalette)) {
        fillPalette <- viridis::plasma(256)
      }

      figure <- figure + geom_point(data = nodes, aes(x = x, y = y, fill = nodes$measure), shape = 21, size = 8, col = "grey") +
        scale_fill_gradientn(limits = fillRange, colours = fillPalette, name = measureName)
    } else {
      figure <- figure + geom_point(data = nodes, aes(x = x, y = y), shape = 21, size = 8, col = "grey",
                                    fill = "sky blue")
    }
  figure <- figure + scale_shape(solid = FALSE) +
    geom_text(data = nodes, aes(x = x, y = y, label = nodes$label), size = label_size) +
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

#' @export
getGraphComponents <- function(obj, threshold = 0.5,
                               minsize = 2, groupNames = NULL) {
  # Finding groups -------------
  library(igraph)
  network <- obj$network
  network[abs(network) < threshold] <- 0
  network[network != 0] <- 1
  graph <- igraph::graph.adjacency(network)
  comp <- igraph::components(graph)
  groups <- which(comp$csize >= minsize)
  groups <- lapply(groups, function(x) names(comp$membership)[comp$membership == x])
  if(is.null(groupNames)) {
    names(groups) <- sapply(groups, function(x) x[1])
  } else if(length(groupNames) != length(groups)) {
    warning("Length of groupNames must equal the number of groups!")
    names(groups) <- sapply(groups, function(x) x[1])
  } else {
    names(groups) <- groupNames
  }

  return(groups)
}
