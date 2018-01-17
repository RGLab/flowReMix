#'@name aggregateModels
#'@title Aggregate flowReMix model fits
#'@description Average flowReMix model fits based on multiple starting points.
#'@param x A \code{list} of flowReMix model fits or a \code{character} vector of \code{rds} files storing model fits.
#'@param verbose \code{logical}. Display a progress bar.
#'@importFrom purrr flatten
#'@importFrom purrr map
#'@export
aggregateModels = function(x, verbose=TRUE, summarizeCoefs = FALSE){
  if(is.list(x)){
      if(!all(unlist(flatten(map(.x = x, .f = function(z)inherits(z,"flowReMix")))))){
        stop("x must be a list of flowReMix model fits.", call. = FALSE)
      }
    TYPE = "list"
  }else if(all(sapply(x,is.character))){
    if(!all(grepl("\\.rds$",x)))
      stop("x must be a vector of rds file names holding flowReMix model fits.", call. = FALSE)
    TYPE = "vector"
  }else{
    stop("x must be a vector of rds file names holing flowReMix model fits or a list of flowReMix model fits.", call. = FALSE)
  }
  output = list()
  coefList = list()
  postList = list()
  levelProbsMatrix = NULL
  models=list()
  if(TYPE=="list"){
    if(verbose){
      pb = txtProgressBar(min = 1, max = length(x),style=3)
    }
    models=x #list of all models
    #average each
    output = x[[1L]]
    coefList = c(coefList,list(output$coefficients))
    postList = c(postList,list(output$posteriors))
    levelProbsMatrix = cbind(levelProbsMatrix,output$levelProbs)
    for(i in seq_along(x)[-1L]){
      if(verbose){
        setTxtProgressBar(pb,i)
      }
      this = x[[i]]
      #we also need the CIs of the coefficients and the posterior probabilities.
      #we'll keep a list of each, then summarize them and save them in a slot in the final object.
      output = .averageModels(output,this,i,x)
      coefList = c(coefList,list(this$coefficients))
      postList = c(postList,list(this$posteriors))
      levelProbsMatrix = cbind(levelProbsMatrix,this$levelProbs)
    }
  }
  if(TYPE == "vector"){
    #load each and average.
    if(verbose){
      pb = txtProgressBar(min = 1, max = length(x),style=3)
    }
    for(i in seq_along(x)){
      this = try(readRDS(x[i]),silent==TRUE)
      models[[i]]=this #keep a list of all models
      if(verbose){
        setTxtProgressBar(pb,i)
      }
      if(inherits(this,"try-error")){
        stop(x[i]," is not a valid rds file", call = FALSE)
      }
      if(!inherits(this,"flowReMix")){
        stop(x[i], " is not a valid flowReMix fit object", call = FALSE)
      }
      #continue merging
      if(i==1){
        output = this
        coefList = c(coefList,list(output$coefficients))
        postList = c(postList,list(output$posteriors))
        levelProbsMatrix = cbind(levelProbsMatrix,output$levelProbs)
      }else{
        output = .averageModels(output,this,i,x)
        coefList = c(coefList,list(this$coefficients))
        postList = c(postList,list(this$posteriors))
        levelProbsMatrix = cbind(levelProbsMatrix,this$levelProbs)
      }
    }
  }
  if(verbose){
    close(pb)
    }
  output$isingCount = round(output$isingCount) #Should be an integer
  #next we compute the CIs

  if(summarizeCoefs) {
    coef_summary = .summarizeCoefs(coefList)
  } else {
    coef_summary <- NULL
  }
  post_summary = .summarizePost(postList)
  rownames(levelProbsMatrix) = colnames(this$posteriors)[-1L]
  levelProbs_summary = .summarizeLevelProbs(levelProbsMatrix)
  output$coef_summary = coef_summary
  output$post_summary = post_summary
  output$coefList = coefList
  output$postList = postList
  output$levelProbs_summary = levelProbs_summary
  output$levelProbsMatrix = levelProbsMatrix
  output$models=models #return a list of all models as part of the averaged fit.
  return(output)
}

#'@importFrom purrr map2
#'@importFrom purrr map2_dbl
#'@importFrom purrr map2_df
#'@importFrom reshape2 melt
#'@importFrom plyr ldply
#'@importFrom tidyr gather
#'@importFrom tidyr spread
.averageModels <- function(output, this, i, x, hasIsing) {
  if(!all.equal(output$posteriors[,1],this$posteriors[,1])|!all.equal(output$data,this$data))
    stop("Models ",x[i], " and ", x[i-1]," are not compatible.")
  #running average of coefficients
  output$coefficients = map2(output$coefficients,this$coefficients,function(x,y)x*(i-1)/i+y*1/i)
  output$covariance = matrix(map2_dbl(output$covariance,this$covariance,function(x,y)x*(i-1)/i+y*1/i),ncol=ncol(output$covariance), dimnames = list(rownames(output$covariance),colnames(output$covariance)))
  output$invCovAvg = matrix(map2_dbl(output$invCovAvg,this$invCovAvg,function(x,y)x*(i-1)/i+y*1/i),ncol=ncol(output$invCovAvg),dimnames = list(rownames(output$invCovAvg),colnames(output$invCovAvg)))
  output$invCovVar = matrix(map2_dbl(output$invCovVar,this$invCovVar,function(x,y)x*(i-1)/i+y*1/i),ncol=ncol(output$invCovVar),dimnames = list(rownames(output$invCovVar),colnames(output$invCovVar)))
  output$randomEffects[,-1L] = as.data.frame(map2_df(output$randomEffects[,-1L],this$randomEffects[,-1L],function(x,y)x*(i-1)/i+y*1/i),check.names=FALSE)
  output$dispersion = map2_dbl(output$dispersion,this$dispersion,function(x,y)x*(i-1)/i+y*1/i)
  output$isingCov = matrix(map2_dbl(output$isingCov,this$isingCov,function(x,y)x*(i-1)/i+y*1/i),ncol=ncol(output$isingCov), dimnames = list(rownames(output$isingCov),colnames(output$isingCov)))
  if(!is.null(output$isingfit)) {
    output$isingfit = matrix(map2_dbl(output$isingfit,this$isingfit,function(x,y)x*(i-1)/i+y*1/i),ncol=ncol(output$isingfit), dimnames = list(rownames(output$isingfit),colnames(output$isingfit)))
  }
  output$posteriors[,-1L] = as.data.frame(map2_df(output$posteriors[,-1L],this$posteriors[,-1L],function(x,y)x*(i-1)/i+y*1/i),check.names=FALSE)
  output$levelProbs = map2_dbl(output$levelProbs,this$levelProbs,function(x,y)x*(i-1)/i+y*1/i)
  output$isingAvg = matrix(map2_dbl(output$isingAvg,this$isingAvg,function(x,y)x*(i-1)/i+y*1/i),ncol=ncol(output$isingAvg), dimnames = list(rownames(output$isingAvg),colnames(output$isingAvg)))
  output$isingVar = matrix(map2_dbl(output$isingVar,this$isingVar,function(x,y)x*(i-1)/i+y*1/i),ncol=ncol(output$isingVar), dimnames = list(rownames(output$isingVar),colnames(output$isingVar)))
  output$isingCount = matrix(map2_dbl(output$isingCount,this$isingCount,function(x,y)x*(i-1)/i+y*1/i),ncol=ncol(output$isingCount), dimnames = list(rownames(output$isingCount),colnames(output$isingCount)))
  if(!is.null(output$isingStability)) {
    output$isingStability$network = matrix(map2_dbl(output$isingStability$network,this$isingStability$network,function(x,y)x*(i-1)/i+y*1/i),ncol=ncol(this$isingCount), dimnames = list(rownames(this$isingStability$network),colnames(this$isingStability$network)))
  }

  #set the class
  if(!inherits(output,"flowReMixAggregate"))
    class(output) = c(class(output),"flowReMixAggregate")
  return(output)
}

.summarizePost = function(postList){
  nr = nrow(postList[[1]])
  nc = ncol(postList[[2]])
  post_without_id = Map(postList,f = function(x){
    x[,-1L]
  })
  post_array = array(unlist(post_without_id),dim = c(nr,nc-1,length(postList)))
  #first dimension is row, second is column, third is model, e.g. post_array[1,3,2] is subject 1, subset 3, model 2
  post_summary = apply(post_array,1:2,function(x).describe(x,df=FALSE))#function(x)c(mean=mean(x),sd=sd(x),n=length(x),quantile(x,c(0.1,0.5,0.9))))
  dimnames(post_summary)[[2]] = rownames(postList[[1]]) # put names on the remaining dimensions, subjects in dim 2 and subsets in dim 3
  dimnames(post_summary)[[3]] = colnames(postList[[1]])[-1L]
  post_summary = melt(aperm(post_summary,perm = c(2,3,1))) #make subjects dim 1, subsets dim 2, and statistics dim 3
  colnames(post_summary)[1:3] = c(colnames(postList[[1]])[1],"subset","statistic")
  post_summary = post_summary%>%spread(statistic,value)
  post_summary
}


.summarizeCoefs <- function(coefList) {
  # coefsummaries = ldply(flatten(coefList)) %>% gather(coef, effect, -.id) %>% group_by(.id, coef) %>%
  #   do({.describe(.$effect)
  #   })
  coefsummaries = Map(
    flatten(coefList),
    f = function(x)
      gather(ldply(x), column, effect, -.id) %>% rename(coef = .id) %>% select(-column)
  ) %>% ldply %>% group_by(.id, coef) %>% do(.describe(.$effect))
  colnames(coefsummaries)[1] = "subset"
  coefsummaries
}

.summarizeLevelProbs = function(levelProbsMatrix){
  levelProbs_summary = do.call(rbind, apply(levelProbsMatrix, 1, .describe))
  rownames(levelProbs_summary) = rownames(levelProbsMatrix)
  levelProbs_summary
}

.describe = function(x,df=TRUE){
  if(df){
    data.frame(
      mean = mean(x),
      sd = sd(x),
      n = length(x),
      t(quantile(x, c(0.1, 0.5, 0.9))),
      check.names = FALSE
    )
  }else{
    c(
      mean = mean(x),
      sd = sd(x),
      n = length(x),
      quantile(x, c(0.1, 0.5, 0.9))
    )
  }
}
