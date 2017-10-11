#'@name aggregateModels
#'@title Aggregate flowReMix model fits
#'@description Average flowReMix model fits based on multiple starting points.
#'@param x A \code{list} of flowReMix model fits or a \code{character} vector of \code{rds} files storing model fits.
#'@param verbose \code{logical}. Display a progress bar.
#'@importFrom purrr flatten
#'@importFrom purrr map
#'@export
aggregateModels = function(x, verbose=TRUE){
  if(is.list(x)){
      if(!all(flatten(map(.x = x, .f = function(z)inherits(z,"flowReMix"))))){
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
  if(TYPE=="list"){
    if(verbose)
      pb = txtProgressBar(min = 1, max = length(x),style=3)
    #average each
    output = x[[1L]]
    for(i in seq_along(x)[-1L]){
      if(verbose)
        setTxtProgressBar(pb,i)
      this = x[[i]]
      output = .averageModels(output,this,i,x)
    }
  }
  if(TYPE == "vector"){
    #load each and average.
    if(verbose)
      pb = txtProgressBar(min = 1, max = length(x),style=3)
    for(i in seq_along(x)){
      this = try(readRDS(x[i]),silent==TRUE)
      if(verbose)
        setTxtProgressBar(pb,i)
      if(inherits(this,"try-error")){
        stop(x[i]," is not a valid rds file", call = FALSE)
      }
      if(!inherits(this,"flowReMix")){
        stop(x[i], " is not a valid flowReMix fit object", call = FALSE)
      }
      #continue work
      if(i==1)
        output = this
      else{
        #check that models are compatible
        output = .averageModels(output,this,i,x)
      }
    }
  }
  if(verbose)
    close(pb)
  output$isingCount = round(output$isingCount) #Should be an integer
  return(output)
}

#'@importFrom purrr map2
#'@importFrom purrr map2_dbl
#'@importFrom purrr map2_df
.averageModels <- function(output,this,i,x) {
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
  output$isingfit = matrix(map2_dbl(output$isingfit,this$isingfit,function(x,y)x*(i-1)/i+y*1/i),ncol=ncol(output$isingfit), dimnames = list(rownames(output$isingfit),colnames(output$isingfit)))
  output$posteriors[,-1L] = as.data.frame(map2_df(output$posteriors[,-1L],this$posteriors[,-1L],function(x,y)x*(i-1)/i+y*1/i),check.names=FALSE)
  output$levelProbs = map2_dbl(output$levelProbs,this$levelProbs,function(x,y)x*(i-1)/i+y*1/i)
  output$isingAvg = matrix(map2_dbl(output$isingAvg,this$isingAvg,function(x,y)x*(i-1)/i+y*1/i),ncol=ncol(output$isingAvg), dimnames = list(rownames(output$isingAvg),colnames(output$isingAvg)))
  output$isingVar = matrix(map2_dbl(output$isingVar,this$isingVar,function(x,y)x*(i-1)/i+y*1/i),ncol=ncol(output$isingVar), dimnames = list(rownames(output$isingVar),colnames(output$isingVar)))
  output$isingCount = matrix(map2_dbl(output$isingCount,this$isingCount,function(x,y)x*(i-1)/i+y*1/i),ncol=ncol(output$isingCount), dimnames = list(rownames(output$isingCount),colnames(output$isingCount)))
  return(output)
}
