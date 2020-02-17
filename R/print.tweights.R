
#' @export
print.tweights <- function(obj) {
  weights=obj$weights
  target=obj$target
  originalTarget=obj$originalTarget
  achievedMean=obj$achievedMean
  Nindependent=obj$Nindependent
  if(is.null(originalTarget)) {
    toprint= t(cbind(achievedMean, target))
    rownames(toprint) =c("Achieved Mean", "Target Mean")
    # colnames(toprint)=colnames(dataset)
  } else {
    toprint= t(cbind(achievedMean, target, originalTarget))
    rownames(toprint) =c("Achieved Mean", "Adjusted Target Mean", "Original Target Mean")
    # colnames(toprint)=colnames(dataset)
  }
  
  
  cat("----------------------------------------------------------------\n")
  cat("Optimization was successful. The weights have a sampleing\ndistribution with means close to the attemted target:\n")
  print(toprint)
  cat("Maximum weight was: ", max(obj$weights),"\n")
  if( Nindependent >0 )
    cat("Data augmented with", Nindependent, "samples with independent variables.",
        "\nThe final weight of these samples was: ", 
        sum(weigths[(length(weights)-Nindependent+1):length(weights)]))
  
  cat("----------------------------------------------------------------\n")
  
}
