#' @name stabilityGraph
#' @title compute a stability graph for the Ising model from a flowReMix fit.
#' @description Compute the stable components of the Ising or Random Effects graphical model through bootstrapping.
#' @param obj \code{flowReMix} model fit.
#' @param type \code{character} "ising" or "randomEffects"
#' @param cv \code{integer} MORE DETAILS
#' @param reps \code{integer} number of replications.
#' @param cpus \code{integer} number of cpus to use
#' @param gamma \code{numeric} MORE DETAILS
#' @param AND \code{logical}. How the Ising edges are aggregated. AND or OR.
#' @param seed random seed. Default 100 if NULL.
#' @import ggplot2
#' @importFrom network network
#' @importFrom GGally ggnet2
#' @importFrom viridis plasma
#' @importFrom igraph graph.adjacency
#' @import doParallel
#' @importFrom parallel splitIndices
#' @importFrom parallel mclapply
#' @export
stabilityGraph <- function(obj, type = c("ising", "randomEffects"),
                           cv = FALSE, reps = 100, cpus = 1,
                           gamma = 0.9, AND = TRUE, seed = NULL,
                           sampleNew = FALSE,keepEach=10) {
  if(cpus == 1) {
    foreach::registerDoSEQ()
  } else {
    registerDoParallel(cores = cpus)
    registerDoRNG()
  }

  if(!is.null(seed)) {
    set.seed(seed)
  }else{
    set.seed(100)
  }

  type <- type[1]
  if(!sampleNew) {
    if(type == "ising") {
      if(is.null(obj$assignmentList)) {
        stop("Posterior samples were not kept, please re-run with `sampleNew = TRUE'.")
      }
      samples <- obj$assignmentList
      names(samples) <- names(obj$randomEffectSamp)
      family <- "binomial"
    } else if(type == "randomEffects") {
      if(is.null(obj$randomEffectSamp)) {
        stop("Posterior samples were not kept, please re-run with `sampleNew = TRUE'.")
      }
      samples <- obj$randomEffectSamp
      family <- "gaussian"
    } else {
      stop("Method not supported!")
    }
  } else {
    mhList <- reConstructMHlist(obj)
    keepEach <- min(obj$control$keepEach, keepEach)
    nsamp <- ceiling(reps * keepEach)
    isingCoefs <- obj$isingAvg
    nSubsets <- length(obj$coefficients)
    intSampSize <- obj$control$intSampSize
    cov <- obj$covariance
    if(!is.null(obj$MHcoef)) {
      MHcoef <- obj$MHcoef
    } else {
      warning("Fit object does not contain a vector of MH coefficients! Using inital (sub-optimal!) values!")
      MHcoef <- rep(obj$control$initMHcoef, nSubsets)
    }
    preCoef <- obj$control$preAssignCoefs[length(obj$control$preAssignCoefs)]
    prior <- obj$control$prior
    dispersion <- obj$dispersion
    invcov <- obj$invCovAvg
    mixed <- is.null(obj$isingAvg)
    # inds <- splitIndices(length(mhList), pmin(length(mhList),cpus))
    inds <- splitIndices(length(mhList), pmin(length(mhList),1))
    mhList <- lapply(inds, function(x) mhList[x])
    mhList = extractFromList(mhList)

    mcEM = obj$control$markovChainEM
    if(!exists("doNotSample")){
      doNotSample = rep(FALSE,nSubsets)
    }
    # if(length(mhList)==1){
      MHresult = CppFlowSstepList_mc_vec(nsubjects = mhList$N, Y = mhList$Y,
                                         N = mhList$TOT, subpopInd = mhList$subpopInd,
                                         clusterassignments = mhList$clusterassignments,
                                         nullEta = mhList$nullEta, altEta = mhList$altEta,
                                         rand = mhList$rand, index = mhList$index, preassign = mhList$preassign,
                                         nsamp = nsamp, nsubsets = nSubsets,
                                         intSampSize = intSampSize, isingCoefs = isingCoefs,
                                         covariance = cov, keepEach = keepEach,
                                         MHcoef = MHcoef, betaDispersion = TRUE,
                                         randomAssignProb = 0,
                                         iterAssignCoef = preCoef, prior = prior,
                                         zeroPosteriorProbs = FALSE,
                                         M = dispersion, invcov = invcov, mixed = mixed,
                                         sampleRandom = (type != "ising"),
                                         doNotSample = doNotSample, markovChainEM = mcEM, cpus=cpus, seed = seed)

    # MHresult =   CppFlowSstepList_mc(mhList[[1]], nsamp = nsamp, nSubsets = nSubsets,
    #                         intSampSize = intSampSize, isingCoefs = isingCoefs,
    #                         covariance = cov, keepEach = keepEach,
    #                         MHcoef = MHcoef, betaDispersion = TRUE,
    #                         randomAssignProb = 0, modelprobs = 0,
    #                         iterAssignCoef = preCoef, prior = prior, zeroPosteriorProbs = FALSE,
    #                         M = dispersion, invcov = invcov, mixed = mixed,
    #                         sampleRandom = (type != "ising"),doNotSample = doNotSample, markovChainEM = mcEM)
    # # }else{
    # MHresult <- foreach(sublist = mhList, .combine = c) %dorng% {
    #   # lapply(sublist, function(subjectData) {
    #   CppFlowSstepList_mc(sublist, nsamp = nsamp, nSubsets = nSubsets,
    #                       intSampSize = intSampSize, isingCoefs = isingCoefs,
    #                       covariance = cov, keepEach = keepEach,
    #                       MHcoef = MHcoef, betaDispersion = TRUE,
    #                       randomAssignProb = 0, modelprobs = 0,
    #                       iterAssignCoef = preCoef, prior = prior, zeroPosteriorProbs = FALSE,
    #                       M = dispersion, invcov = invcov, mixed = mixed,
    #                       sampleRandom = (type != "ising"),doNotSample = doNotSample, markovChainEM = mcEM)
    #      }
    # }
    if(type == "ising") {
      # samples <- lapply(MHresult, function(x) x$assign)
      samples = MHresult$assign
      family = "binomial"
    } else if(type == "randomEffects") {
      # samples <- lapply(MHresult, function(x) x$rand)
      samples = MHresult$rand
      family <- "gaussian"
    }
  }

  if(!sampleNew) {
    #what if the samples don't have names?
    #TODO fix
    names(samples) <- sapply(names(samples), function(x) strsplit(x, "%%%")[[1]][[1]])
    samples <- lapply(unique(names(samples)), function(x) {
      do.call("rbind", samples[names(samples) == x])
    })
  }
  subsets <- names(obj$coefficients)
  nsubsets <- ncol(samples[[1]])

  if(sampleNew) {
    #no need for this since we are already in the right format using mc_vec function
    # samples <- lapply(1:reps, function(i) t(sapply(samples, function(samp) samp[i, , drop = FALSE])))
    #below is 20 times faster.
    #taking slices along the third dimension gives contiguous memory chunks..should do an aperm so that reps is the third dim.
    # samples = array(unlist(samples,use.names=FALSE),dim = c(nrow(samples[[1]]),ncol(samples[[1]]),length(samples)))
    # samples = aperm(samples,c(3,2,1))
  } else {
    #TODO fix
    # samples <- lapply(1:reps, function(i) t(sapply(samples, function(samp) samp[sample.int(nrow(samp), 1), , drop = FALSE])))
    samples2 = array(unlist(samples,use.names=FALSE),c(nrow(samples[[1]]),ncol(samples[[1]]),length(samples))) #make it a multidimensional array
    samples = array(apply(samples2,3,function(x)x[sample.int(nrow(x)),]),dim=c(nrow(samples[[1]]),ncol(samples[[1]]),length(samples))) #and permute the replicates per subject
    # samples = aperm(samples,c(3,2,1)) #TODO check that it should be correct so we to t(mat) below..
  }

  if(!sampleNew){
    inds = splitIndices(reps,pmin(reps,cpus))
    for(i in 1:length(inds)){
      inds[[i]]=rep(i,length(inds[[i]]))
    }
    inds = split(sample(dim(samples)[1],reps,replace=TRUE),inds)
  }else{
    inds = splitIndices(reps, pmin(reps,cpus))
  }
  #splitting third dimension below.. is along replicates.. so we change that to split by the first now that we've permuted
  samples <- lapply(inds, function(x) samples[x,,,drop=FALSE])
  set.seed(seed)
  countCovar <- foreach(matrices = samples , .combine = "+") %dorng% {
    countCovar <- apply(matrices,1, function(mat) #by first dim - replicates
      raIsing(t(mat), AND = AND, gamma = gamma, family = family,
              method = "sparse", cv = cv, parallel = FALSE) != 0)
    countCovar = matrix(rowSums(countCovar),ncol=dim(matrices)[2])
    return(countCovar)
  }
  # countCovar <- Reduce(cluster_res, f = "+")
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

# Function for reconstructing MHlist for sampling new posteriors --------------
reConstructMHlist <- function(obj) {
  mixed <- is.null(obj$isingfit)
  dat <- buildFlowFrame(obj$call, obj$data)
  dat <- dat[[1]]
  formulas <- editFlowFormulas(obj$call, mixed)
  initFormula <- formulas$initFormula
  dat$iteration <- 1
  bypop <- by(dat, INDICES = dat$sub.population, FUN = computeFlowEta, obj$coefficients, 1, initFormula, mixed)
  names(bypop) <- names(obj$coefficients)
  databyid <- do.call("rbind", bypop)
  rm(bypop)
  databyid <- with(databyid, databyid[order(sub.population, id, decreasing = FALSE), ])
  databyid <- by(databyid, databyid$id, function(x) x)
  preAssignment <- do.call("rbind", by(dat, dat$id, autoPreAssign))
  preAssignment <- preAssignment[order(preAssignment$id, preAssignment$subset), ]
  preAssignment <- by(preAssignment, preAssignment$id, function(x) x)

  forcols <- databyid[[1]]
  keepcols <- which(names(forcols) %in% c("y", "N", "subpopInd", "nullEta", "altEta"))
  rm(forcols)

  estimatedRandomEffects <- obj$randomEffects
  if(class(estimatedRandomEffects)=="data.frame"|ncol(estimatedRandomEffects)>length(nlevels(factor(dat$sub.population)))){
    estimatedRandomEffects = as.matrix(estimatedRandomEffects[,-1L])
  }
  nSubjects <- nrow(obj$posteriors)
  listForMH <- lapply(1:nSubjects, function(i, keepcols) list(dat = databyid[[i]][, keepcols],
                                                              pre = preAssignment[[i]],
                                                              rand = estimatedRandomEffects[i, ],
                                                              index = i), keepcols)
  names(listForMH) <- names(databyid)
  return(listForMH)
}

plotRawGraph <- function(obj, graph = c("ising"), threshold = 0.5, plotAll = FALSE,
                         fill = NULL, fillName = NULL,
                         fillRange = NULL, fillPalette = NULL,
                         title = TRUE, normalize = FALSE,
                         count = TRUE, label_size = 1.8,seed=100, ...) {
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
  rownames(isingZ) <- names(obj$coefficients)
  colnames(isingZ) <- names(obj$coefficients)
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
#' @import ggplot2
#' @export
plot.flowReMix_stability <- function(x, threshold = 0.5, nEdges = NULL, plotAll = FALSE,
                                     fill = NULL, fillName = NULL, fillRange = NULL,
                                     fillPalette = NULL, title = TRUE,
                                     layout = "fruchtermanreingold",
                                     layout.par = NULL, label_size = 1.8, seed = 100,
                                     ...){
  set.seed(seed)
  measure <- fill
  props <- x$network
  if(is.null(nEdges)) {
    props[abs(props) < threshold] <- 0
  } else {
    propVals <- sort(unique(as.vector(abs(props))), decreasing = TRUE)
    propCounts <- sapply(propVals, function(x) sum(props == x) / 2)
    threshold <- propVals[min(which(cumsum(propCounts) >= nEdges))]
    props[abs(props) < threshold] <- 0
  }
  network <- props
  if(!plotAll) {
    keep <- apply(network, 1, function(x) any(abs(x) >= threshold))
  } else {
    keep <- rep(TRUE, ncol(network))
  }
  network <- network[keep, keep]
  net <- network::network(props)
  subsets <- colnames(props)
  nodes <- GGally::ggnet2(network, label = subsets[keep], mode = layout,
                          layout.par = layout.par)$data
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

  if(is.null(x$raw)) {
    alphaTitle <- "Edge Propensity"
  } else {
    alphaTitle <- "Edge Quantile"
  }

  requireNamespace("ggplot2")
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
  } else if(x$type == "ising" & title) {
    figure <- figure + ggtitle("Estimated Dependence Structure for Ising Model")
  } else if(x$type == "randomEffects" & title) {
    figure <- figure + ggtitle("Estimated Dependence Structure for Random Effects")
  }

  return(figure)
}

#' @importFrom igraph graph.adjacency
#' @importFrom igraph components
getGraphComponents <- function(obj, threshold = 0.5,
                               minsize = 2, groupNames = NULL) {
  # Finding groups -------------
  network <- obj$network
  network[abs(network) < threshold] <- 0
  network[network != 0] <- 1
  graph <- graph.adjacency(network)
  comp <- components(graph)
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
