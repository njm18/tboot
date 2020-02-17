#' @title Function tboot
#' @description Bootstrap \code{nrow} rows of \code{dataset} using
#' the given row-level weights.
#' @seealso \code{\link{tweights}}
#' @export
#' @param nrow number of rows in the new bootstrapped dataset.
#' @param weights an object of class 'tweights' output from the 'tweights' function.
#' @param dataset Data frame or matrix to bootstrap. Rows of the dataset must be in the 
#' same order as was used for the 'tweights' call. However the dataset may include
#' additional columns not included in the 'tweights' calll.


#' @param fillMissingAug fill in missing augmentation with primary weights resampling.
tboot <- function(nrow,
                  weights,
                  dataset=weights$dataset,
                  fillMissingAug=TRUE) {
  if(missing(nrow))
    stop("'nrow' is missing")
  if(missing(weights))
    stop("'weights' is missing")
  
  nweights <- length(weights$weights)
  index <- sample.int(
    n = nweights,
    size = nrow,
    prob = weights$weights,
    replace = TRUE
  )
  
  if(!("tweights" %in% class(weights)))
    stop("'weights' must be an object of class 'tweights' from the 'tweights' function.")
  
  Nindependent=weights$Nindependent
  if(is.null(weights$Nindependent))
    Nindependent=0
  if(Nindependent==0) {
    if (nweights != nrow(dataset)){
      stop("length of weights must be nrow(dataset).")
    }
    return(dataset[index, ,drop=FALSE])
  } else {
    
    #need to deal with augmentation
    if(nweights != (Nindependent+nrow(dataset)))
      stop("length of weights must be nrow(dataset)+Nindependent.")
    

    if(is.null(weights$augmentWeights))
      stop("Attributes of weights not set correctly for 'augmentWeights.'")
    if( !any(class(weights$augmentWeights) =="list") )
      stop("'augmentWeights' must be a 'list.'")
    if( is.null(names(weights$augmentWeights)) )
      stop("'augmentWeights' must be a named 'list.'")
    
    #Fill in any missing in case an unconstrained variable was added
    missingAug=names(dataset)[!(names(dataset) %in% names(weights$augmentWeights))]
    if(length(missingAug)>0) {
      if(!fillMissingAug) {
        stop("Missing 'augmentWeights.' Consider setting fillMissingAug=TRUE.")
      } 
      wtmp=weights[1:(length(weights)-Nindependent)]
      wtmp=wtmp/sum(wtmp)
      for(nm in missingAug) {
        weights$augmentWeights[[nm]]=wtmp
      }
    }
    
    #augmented samples start as NA
    index[index > nrow(dataset)]=NA
    
    #get non-augmented
    ret=dataset[index, ,drop=FALSE] 
    #cat(which(is.na(ret[,1])))
    
    #get augmented
    augIndex=which(is.na(index))
    naug=length(augIndex)
    aug=do.call(data.frame, lapply(colnames(dataset), 
               function(nm) {
                 w=weights$augmentWeights[[nm]]
                 cat(nm,"\n")
                 cat(length(w),"\n")
                 cat(nrow(dataset),"\n")                 
                 if(length(w)!=nrow(dataset))
                   stop("'augmentWeights' weights not set to correct length.")
                 augindex=sample.int(n = length(w),
                                     size = naug, prob = w, 
                                     replace = TRUE)
                 dataset[augindex,nm]
               }))
    
    ret[augIndex,]=aug
                          
    return(ret)
  }
} 

