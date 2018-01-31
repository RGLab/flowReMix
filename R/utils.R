#'@importFrom stats model.frame model.matrix na.omit na.pass quantile sd
#'@importFrom utils getFromNamespace
NULL

.isFlowRemix = function(x){
  if(!inherits(x,"flowReMix")){
    stop("x must be a flowReMix object.",call. = FALSE)
  }
}

#' getPosteriors
#'
#' @param x \code{flowReMix} model fit
#'
#' @return \code{data.frame} of posterior probabilities for each subject and subset in the model.
#' @export
#'
#' @examples
#' data(fit505)
#' getPosteriors(fit505)
#'
getPosteriors = function(x){
  .isFlowRemix(x)
  return(x$posteriors)
}

#' getSubsets
#'
#' @param x \code{flowReMix} object
#'
#' @return \code{vector} of subset names fitted in the model.
#' @export
#'
#' @examples
#' data(fit505)
#' getSubsets(fit505)
#'
getSubsets = function(x){
  .isFlowRemix(x)
  return(colnames(x$posteriors)[-1L])
}

#' getIsing
#'
#' @param x \code{flowReMix} model fit.
#'
#' @return \code{isingStability} slot from the model.
#' @export
#'
#' @examples
#' data(fit505)
#' getIsing(fit505)
getIsing = function(x){
  .isFlowRemix(x)
  return(x$isingStability)
}


#' degreeFromStringFun
#' A function to parse the degree of a subset from its name.
#' @details Assumes that a cell subset name follows a specific pattern where
#' the positive functions are simply listed in the name. For example, "8+/A+B+C+" would be  a degree 3
#' subset of CD8+ cells.
#' @param x \code{character} cell subset name from \code{getSubsets}
#'
#' @return \code{numeric} vector of degree of functionality.
#' @export
#' @seealso \code{\link{flowReMixPFS}} \code{\link{weightForPFS}}
#' @examples
#' data(fit505)
#' degreeFromStringFun(getSubsets(fit505))
#'
degreeFromStringFun = function(x,split="\\+"){
  degree = unlist(
    Map(
      strsplit(x = gsub(".*?/", "", x), split = split),
      f=length),
    use.names = FALSE)
  names(degree) = x
  degree
}

#' weightForPFS
#' Compute the weighting for the polyfunctionality score for each type of subset.
#' @details Computes the weight for the polyfunctionality score calculation. The weight is \code{choose(M,k)} where
#' k is the degree of a subset and M is the \code{max(degree)} or the number of markers.
#' @param x \code{flowReMix} fit object
#' @param M \code{numeric} the max possible functionality. Uses \code{max(degree)} if not provided.
#' @param parsefun \code{function} to parse the degree from a cell subset name. Defaults to \code{\link{degreeFromStringFun}} in the package.
#' @param split \code{character} used to split the functionality and extract the degree.
#'
#' @return \code{numeric} vector of weights.
#' @seealso \code{\link{flowReMixPFS}} \code{\link{degreeFromStringFun}}
#' @export
#'
#' @examples
#' data(fit505)
#' weightForPFS(fit505)
weightForPFS = function(x, M = NULL, parsefun = degreeFromStringFun, split = "\\+"){
  if("character"%in%class(x)){
    degree =  parsefun(x, split = split)
  }else{
    .isFlowRemix(x)
    degree = parsefun(getSubsets(x),split = split)
  }
  if(is.null(match.call()$M))
     M = max(degree)
  return(degree/choose(M,degree))
}


#' perSubsetPFS
#' Compute the PFS (polyfunctionality score) per subject and stimulation group
#' @param x \code{flowReMix} model object
#' @param M \code{numeric} the max possible cell subset degree, the number of functions measured.
#' @param stimVar \code{name} Unquoted  name of the stimulation variable in the data. e.g. stim
#' @param parentVar \code{name} Unquoted name of the parent cell population variable in the data. e.g. parent
#' @param outcomeVar \code{name} Unquoted name of the outcome variable in the data. If provided, the scores will be merged with the data using the \code{subject_id}
#' @param ... additional arguments passed to weightForPFS.
#' @details
#' Requires that the data table has a variable for stimulation and for cell population parent.
#'
#' The user can pass in a function to the argument  \code{parser=}. This function parses a cell subset name string  and
#' returns the degree of functionality for the cell subset. The default is \code{\link{degreeFromStringFun}}
#' @note Uses rlang quosures to pass variable information. PFS are NOT normalized to the total number of possible subsets (n*(n+1)/2), like in COMPASS.
#' @seealso \code{\link{degreeFromStringFun}} \code{\link{weightForPFS}}
#' @return \code{data.frame} of polyfunctionality scores for each cell subset and subject, weighted appropriately.
#' @export
#'
#' @examples
#' data(fit505)
#' flowReMixPFS(fit505,M=5, stimVar = stimGroup, parentVar = parent)
#' flowReMixPFS(fit505,M=5,stimVar=stimGroup,parentVar=parent,outcomeVar=infection)
#'\dontrun{
#' flowReMixPFS(rv144_aggregate,M=6,stimVar=stim,parentVar=parent,split=",")
#' }
flowReMixPFS = function(x,M,stimVar = NULL, parentVar = NULL, outcomeVar = NULL, ...){
  mc  =  match.call()
  if(is.null(mc$stimVar)|is.null(mc$parentVar)){
    stop("stimVar and parentVar must be provided")
  }
  stimVar = enquo(stimVar)
  parentVar = enquo(parentVar)
  post = getPosteriors(x) %>% gather(pop,posterior,-1) %>% mutate(weight=weightForPFS(pop,M=M, ...))
  if(length(unique(post$weight))==1){
     if(class("split")=="character"){
	split="\\+"
     }
     warning(paste0("All PFS weights are the same! Splitting functionality by ",split," : verify that this is correct!"),call. = FALSE)
  }
  post$score = post$weight*post$posterior
  subjid = x$subject_id
  subjid = enquo(subjid)
    post = x$data %>%
      select(!!parentVar,pop=subset,!!stimVar) %>%
      mutate(pop=as.character(pop)) %>%
      unique() %>%
      inner_join(post) %>%
      group_by(!!parentVar,!!stimVar, !!subjid) %>%
      summarize(PFS=mean(score))
    if(!is.null(mc$outcomeVar)){
      outcomeVar = enquo(outcomeVar)
      outcome = x$data%>%select(!!subjid,!!outcomeVar)%>%unique()
      post = inner_join(outcome,post)
    }

  return(post)
}


#' plotIsingGraph
#'
#' Plot an ising graph.
#' @param x \code{flowReMix} object.
#'
#' @return \code{ggplot} ggraph object
#' @export
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom tidygraph as_tbl_graph
#' @importFrom tidygraph group_components
#' @importFrom tidygraph activate
#' @importFrom rlang enquo
#' @examples
#' data(fit505)
#' plotIsingGraph(fit505,weight=0.9) +
#' guides(shape=guide_legend(nrow=1),
#'   size=guide_legend(nrow=1),
#'   color=guide_legend(nrow=4))
plotIsingGraph = function(x, weight=0.6,layout="kk"){
  weight = enquo(weight)
  .isFlowRemix(x)
  ising = getIsing(x)
  ising_graph = graph_from_adjacency_matrix(getIsing(x)$network, weighted = TRUE,mode = "undirected",diag = FALSE)
  ising_tidy = as_tbl_graph(ising_graph)
  ising_tidy %>% activate(edges) %>% filter(weight > !!weight) %>% activate(nodes) %>%
    filter(name %in% name[.E()$from] |
             name %in% name[.E()$to]) %>%
    mutate(community = as.factor(group_components())) %>%
    mutate(
      degree = (degreeFromStringFun(name)),
      stim = do.call(rbind, strsplit(name, "/"))[, 1],
      parent = do.call(rbind, strsplit(name, "/"))[, 2]
    ) %>%
    ggraph(layout = layout) +
    geom_edge_link() +
    geom_node_point(aes(
      color = community,
      size = degree,
      shape = stim
    )) + theme_graph() + geom_node_text(aes(label=parent),repel=TRUE)
}

#' plotChordDiagram
#' Generate a chord diagram from a flowReMix fit.
#' @details Parses the \code{subsets} labels into factors that can be used to annotate tracks on a chord diagram.
#' Connections correspond to edges in the ising stability graph with a minumum \code{threshold} propensity.
#' A label such as "env/8+/IL2+" with an argument of \code{separator="/"} and \code{varnames=c("stim","parent","functions")} would
#' cause the cell subset label to be parsed into three variables corresponding to the stimulation, the
#' parent cell population, and the cytokine function.
#' @param x \code{flowReMix} fit
#' @param threshold \code{numeric} threshold for the minimum edge weight in the graph (exclusive).
#' @param varnames \code{character} variable names to define by splitting subsets. e.g. \code{c("parent","stim","function")}
#' @param separator \code{character} separator to split subset string into variables.
#' @param function.separator \code{character} separator to split multifunction string  into individual functions.
#' @param track.order \code{character} vector of the same length as \code{varnames} specifying the order of the variables in the tracks.
#' @param sort.order \code{character} vector of the same length as \code{varnames} specifying the order of the variables in the tracks.
#' @param inner.track \code{character} specifying which variable in \code{varnames} is on the inner track.
#' @param track.height \code{numeric} vector of track heights, length 1 or \code{length(varnames)}
#' @param track.margin \code{numeric} vector of length 2 specifying the track margin.
#' @param functions.name \code{character} length 1, name of the column holding the functions.
#' @param wrap.at \code{character} regex, string position to wrap labels for all tracks.
#' @param label.cex \code{numeric} scaling factors for the labels in each track.Should be same length as \code{varnames}.
#' @param facing \code{character} passed to circos.text. \code{facing=c("inside", "outside", "reverse.clockwise", "clockwise",
#' "downward", "bending", "bending.inside", "bending.outside")}
#' @param parsefun \code{function} that parses multiple function labels to single function matrix.
#' @param functionality.filter \code{numeric} lower threshold for functionality degree. Exclusive.
#' @param auc.filter \code{numeric} lower threshold for AUC prediction. Exclusive.
#' @param auc.outcome \code{name} unquoted variable name used for prediction to calculate the AUC.
#' @param response.filter \code{numeric} minimum cell subset response probability to include.
#' @return NULL
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom tidygraph as_tbl_graph
#' @importFrom tidygraph group_components
#' @importFrom tidygraph .E
#' @importFrom tidyr separate
#' @importFrom tidyr unnest
#' @importFrom tidyr spread
#' @importFrom rlang sym
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom scales brewer_pal
#' @import circlize
#' @export
#'
#' @examples
#' data(fit505)
#' #also defined in the package namespace
#' parser = function(x,separator,functions){
#'   functions=enquo(functions)
#'    dat = x%>% mutate(func = strsplit(!!functions, split = separator)) %>%
#'   unnest()%>%tidyr::spread(func,func)
#'  dat[is.na(dat)]=""
#'  dat
#' }
#' plotChordDiagram(fit505,
#'    threshold = 0.8,
#'    varnames=c("stim","parent","functions"),
#'    track.order = c("stim","parent","functions"),
#'    sort.order=c("parent","stim","degree"),
#'    functions.name="functions",
#'    track.height = 2.5,
#'    track.margin = c(0,0),
#'    label.cex=0.6,
#'    wrap.at = " ",
#'    functionality.filter=1,
#'    auc.filter = 0,
#'    auc.outcome = hiv,
#'    response.filter = 0.05,
#'    facing=c("inside","inside","inside"),parsefun=parser)
#'    \dontrun{
#'    plotChordDiagram(
#'      rv144_aggregate,
#'      threshold = 0.8,
#'      varnames = c("functions"),
#'      separator = " ",
#'      function.separator = ",",
#'      track.order = "functions",
#'      sort.order = c("degree", "functions"),
#'      inner.track = "functions",
#'      parsefun = parser,
#'      track.height = 2,
#'      track.margin = c(0, 0),
#'      functions.name = c("functions"),
#'      facing = "inside",
#'      wrap.at = " ",
#'      label.cex = 0.6,
#'      functionality.filter=1
#' )
#'    }
plotChordDiagram  = function(x,threshold = 0.6,
                             varnames = c("stim","parent","functions"),
                             separator = "/",
                             function.separator = "\\+",
                             track.order = c("functions","parent","stim"),
                             sort.order = c("functions","parent","stim"),
                             inner.track = "parent",
                             track.height = 2,
                             track.margin = c(0,0),
                             functions.name="functions",
                             wrap.at = c("\\+","",""),
                             label.cex = c(0.5,0.5,0.5),
                             facing = "clockwise",
                             parsefun=NULL,
                             functionality.filter=NULL,
                             auc.filter = NULL,
                             auc.outcome=NULL,
                             response.filter=NULL){

  nw = flowReMix::getIsing(x)$network
  aucfilt = rep(TRUE,ncol(nw))
  if(!is.null(match.call()$auc.filter)){
    auc.outcome = enquo(auc.outcome)
    auc.outcome = sym(quo_name(auc.outcome))
    auctable = eval(as.call(list(getFromNamespace("summary.flowReMix",ns = "flowReMix"),object=x,target=auc.outcome)))%>%
      arrange(-auc)%>%filter(auc > auc.filter)
    aucfilt = colnames(nw)%in%auctable$subset
  }

  degfilt = rep(TRUE,ncol(nw))
  if(!is.null(match.call()$functionality.filter)){
    deg = flowReMix::degreeFromStringFun(x = colnames(nw),split = function.separator)
    degfilt = deg>functionality.filter
  }

  respfilt = rep(TRUE,ncol(nw))
  if(!is.null(match.call()$response.filter)){
    lp = x$levelProbs
    names(lp) = colnames(nw)
    respfilt = lp>response.filter
  }
  nw = nw[respfilt&aucfilt&degfilt,respfilt&aucfilt&degfilt]

  if(!functions.name%in%varnames){
    stop("functions.name not in varnames",call. = FALSE)
  }

  gr = igraph::graph_from_adjacency_matrix(nw,weighted = TRUE,mode = "undirected",diag = FALSE)
  gr = as_tbl_graph(gr)
  nodes = gr %>% activate(nodes) %>% as.data.frame()
  edges = gr %>% activate(edges) %>% as.data.frame()
  edgelist = inner_join(nodes %>%
                          mutate(from = as.integer(rownames(.))), edges) %>%
    mutate(from = NULL) %>% rename(from = name) %>%
    inner_join(nodes %>% mutate(from = as.integer(rownames(.))) %>%
                 rename(to = from)) %>%
    mutate(to = NULL) %>% rename(to = name) %>%
    select(from,to,weight)

      edgelist = edgelist %>% separate(from,
                                   paste0(varnames, "_from"),
                                   sep = separator, remove = FALSE) %>%
        separate(to,
                        paste0(varnames, "_to"), sep = separator, remove = FALSE)
  fromvars = paste0(varnames,"_from")
  tovars = paste0(varnames,"_to")
  edgelist[is.na(edgelist)]=""
    #Map to interaction variables
    edgelist = edgelist %>% mutate(sg_from = interaction(sep = "\n", lapply(fromvars, function(x)
      eval(as.name(
        x
      )))),
      sg_to = interaction(sep = "\n", lapply(tovars, function(x)
        eval(as.name(
          x
        ))))) %>%
      filter(weight > threshold)

  node_clusters = gr %>%
    activate(edges) %>%
    filter(.E()$weight > threshold) %>%
    activate(nodes) %>%
    filter(name %in% name[.E()$from] |name %in% name[.E()$to]) %>%
    mutate(cluster = group_components()) %>% activate(nodes) %>% as.data.frame()


  circos.clear()
  indices = gsub(separator, "\n", node_clusters$name)
  nclust = length(unique(node_clusters$cluster))
  #this needs to be defined in the input
  attribs = do.call(rbind, strsplit(indices, "\n"))
  colnames(attribs) = varnames
  attribs = cbind(name=indices,attribs)
  qfunctions.name=quo(eval(as.name(functions.name)))
  attribs = as.data.frame(attribs,stringsAsFactors = FALSE)%>%parsefun(functions=!!qfunctions.name,separator=function.separator)
  rest = attribs[,setdiff(colnames(attribs),c("name",functions.name))]

  dat = (edgelist)[, c("sg_from", "sg_to", "weight")]%>%mutate(sg_from=as.character(sg_from),sg_to=as.character(sg_to))
  dat = node_clusters %>% mutate(name = gsub(separator, "\n", name)) %>% rename(sg_from = name) %>% inner_join(dat)
  clust = dat$cluster


  tosortby = cbind(attribs,degree=degreeFromStringFun(attribs$name,split=function.separator))[,c(sort.order)]
  o = eval(as.call(c(order,as.list(tosortby))))
  .pad = function(x,l){
    if(length(x)==1){
      x = rep(x,l)
    }else if(length(x)!=l){
      x = rep(x,length.out=l)
    }
    x
  }
  #Make a list of tracks
  .makeTracks = function(ntracks = NULL, track.height = NULL,track.margin = NULL){
    if(ntracks != length(track.height)){
      stop("track.height must be the same length as the varnames: ",ntracks,"\n",call. = FALSE)
    }
    if(length(track.margin)!=2){
      stop("track.margin must be the length 2: ",ntracks,"\n",call. = FALSE)
    }

    tracks=list()
    for(i in seq.int(ntracks)){
      tracks[[i]] = list(track.height = uh(track.height[i], "mm"),track.margin = uh(track.margin, "mm"))
    }
    tracks
  }
  mytracks = .makeTracks(ntracks=ncol(rest),track.height = .pad(track.height,ncol(rest)),track.margin = track.margin)


  #colors for the inner grid and connections.
  inner.grid.colors = colorRampPalette(brewer_pal(type = "qual", palette = 3)(12))(nlevels(factor(attribs[,inner.track])))[factor(attribs[,inner.track])]
  connection_colors = colorRampPalette(brewer_pal(type = "div", palette =
                                6)(8))(nclust)[factor(clust)]
  names(inner.grid.colors) = names(indices)
  facing = .pad(facing,(ncol(attribs)))
  track.height = .pad(track.height,(ncol(attribs)))
  wrap.at = .pad(wrap.at,(ncol(attribs)))
  label.cex = .pad(label.cex, ncol(attribs))
  chordDiagram(dat %>% select(sg_from, sg_to, weight),
               directional = 0,
               order = indices[o],
               direction.type = c("diffHeight"),
               symmetric = TRUE, self.link = 0, annotationTrack = NULL,
               preAllocateTracks = mytracks,
               grid.col = inner.grid.colors[o],
               col = connection_colors,
               link.sort = TRUE)

  #colors for the remianing variables
  track.order = setdiff(track.order,functions.name)
  rest = rest[,c(track.order,setdiff(colnames(rest),track.order))]
  sapply(1:ncol(rest),function(columnidx){
    # browser()
    this.track = (split(get.all.sector.index(), rest[o,columnidx]))
    this.tracklength = length(this.track)
    #decide if we interpolate the color palette
    mxc = (subset(brewer.pal.info,category=="qual")[(columnidx%%8+1),,drop=FALSE]$maxcolors)[[1]]
    track.colors = colorRampPalette(brewer_pal(type = "qual", palette = (columnidx%%8+1))(max(1,mxc-1)))(this.tracklength)
    names(track.colors)=names(this.track)
    track.colors[names(track.colors)%in%""]="#FFFFFF"
    mapply(
      this.track,
      FUN = highlight.sector,
      track.index = columnidx,
      col = track.colors,
      text = gsub(paste0("(",wrap.at[columnidx],")"),"\\1\n",names(this.track)),
      cex = label.cex[columnidx],
      text.col = "black",
      niceFacing = TRUE,
      facing = facing[columnidx]
    )
  })
  invisible(NULL)
}


#' Cell subset label parser
#' Parses cell subset labels into matrix of functions.
#' @param x \code{data.frame0} input with a column corresponding to the cell subset function labels
#' @param separator \code{string} to split the label on
#' @param functions \code{name} bare, unquoted name of the column holding the function label
#' @description Best explained via the example.
#' @seealso \code{\link{plotChordDiagram}}
#' @importFrom tidyr unnest
#' @importFrom tidyr spread
#' @return matrix of individual functions.
#' @export
#'
#' @examples
#' dat = data.frame(name=c("IFNg+TNFa+IL2+","IFNg+","IL2+","CD154+TNF+"),stringsAsFactors=FALSE)
#' parser(dat,separator = "\\+",functions=name)
parser = function(x,separator,functions){
    functions=enquo(functions)
     dat = x%>% mutate(func = strsplit(!!functions, split = separator)) %>%
   unnest()%>%spread(func,func)
   dat[is.na(dat)]=""
   dat
}

completeVec <- function(vec, maxlength, impute = -1) {
  if(length(vec) == maxlength) {
    return(vec)
  } else {
    return(c(vec, rep(impute, maxlength - length(vec))))
  }
}

extractFromList = function(mhlist){
  N = length(mhlist[[1]])
  J = nrow(mhlist[[1]][[1]]$pre)
  Y <- map(map(mhlist[[1]],"dat"),"y")
  maxlength <- max(sapply(Y, length))
  Y <- matrix(unlist(lapply(Y, completeVec, maxlength)), ncol = N)
  TOT = matrix(unlist(lapply(map(map(mhlist[[1]],"dat"),"N"), completeVec, maxlength)), ncol = N)
  subpopInd = matrix(unlist(lapply(map(map(mhlist[[1]],"dat"),"subpopInd"), completeVec, maxlength)), ncol = N)
  nullEta = matrix(unlist(lapply(map(map(mhlist[[1]],"dat"),"nullEta"), completeVec, maxlength)), ncol = N)
  altEta = matrix(unlist(lapply(map(map(mhlist[[1]],"dat"),"altEta"), completeVec, maxlength)), ncol = N)
  rand = matrix(unlist(map(mhlist[[1]],"rand")),ncol=N)
  index = unlist(map(mhlist[[1]],"index"))
  preassign = matrix(unlist(map(map(mhlist[[1]], "pre"),"assign")), ncol=N)
  if(length(mhlist[[1]][[1]]$assign) != 0){
    clusterassignments = matrix(unlist(map(mhlist[[1]], "assign")), ncol=N);
  }else{
    clusterassignments = matrix(0, nrow = J, ncol = N)
  }
  return(list(N=N, J=J,Y=Y,TOT=TOT,subpopInd = subpopInd,nullEta=nullEta,altEta=altEta,rand=rand,index=index,preassign=preassign,clusterassignments=clusterassignments))
}

#'@name zeroPosteriorProbs
#'@title zeroPosteriorProbs
#'description return the posterior probabilities with preassignments set to zero.
#'@param modelfit A flowReMix model fit
#'@return matrix of posterior probabilities
#'@export
zeroPosteriorProbs = function(modelfit){
  pre = modelfit$preAssignment
  posteriors = reshape2::melt(getPosteriors(modelfit),value.name="posterior",variable.name="subset",id=quo_name(modelfit$subject_id))

  colnames(posteriors)[1]="id"
  pre$id  = unlist(map(.x = strsplit(as.character(pre$id),"%%%"),.f=function(x)x[1]))
  posteriors$id = factor(posteriors$id)
  pre$id = factor(pre$id)
  zpost = pre%>%unique%>%inner_join(posteriors)%>%mutate(zpost = ifelse(assign==0,0,posterior))
  zpost = zpost%>%select(-assign,-posterior)%>%reshape2::dcast(id~subset,value.var="zpost")
  colnames(zpost)[1]=as.character(modelfit$subject_id)
  return(zpost)
}

#' Translate COMPASS markere names to FlowReMix format
#'
#' @param x \code{matrix} with column names in COMPASS format
#'
#' @return \code{matrix} with column names in FlowReMix format.
#' @export
#'
#' @examples
#' \dontrun{
#' #x is a matrix with columns named by COMPASS.
#' translateCOMPASStoFlowReMixNames(x)
#' }
translateCOMPASStoflowReMixNames = function(x){
   frmmn = map(.x=strsplit(gsub("!.*?$","",gsub(perl=TRUE,"!.*?[&$]","",colnames(x))),"&"),.f=function(x)paste0(paste(x,collapse="+"),"+"))
   colnames(x) = frmmn
   return(x)
 }

#' Create a table suitable for FlowReMix modeling from a list of COMPASS container objects
#'
#' @param x \code{list} of COMPASSContainer objects
#'
#' @return \code{data.frame} suitable for input to flowReMix.
#' @export
#'
#' @examples
#' \dontrun{
#' #x is a list of COMPASSContainer
#' FlowReMixFromCOMPASSContainers(x)
#'
#' }
FlowReMixFromCOMPASSContainers = function(x){
  counts = map(.x = x,.f=CellCounts)
   counts = map(.x=counts,translateCOMPASStoflowReMixNames)
   counts = map(counts,reshape2::melt)
   counts = plyr::ldply(counts)
   colnames(counts) = c("parent",x[[1]]$sample_id,"population","count")
   counts = inner_join(counts,unique(do.call(rbind,map(x,function(x)x$meta))))
   totals = unique(ldply(map(map(x,function(x)function(y=as.data.frame(x$counts)){data.frame(sample_id=rownames(y),count=y)}),function(x)x())))
   colnames(totals) = c("parent",x[[1]]$sample_id,"parentcount")
   counts = inner_join(totals,counts)
   counts
 }

