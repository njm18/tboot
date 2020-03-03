#' @method print tweights    
#' @export
print.tweights <- function(x, ...) {
  weights=x$weights
  target=x$target
  originalTarget=x$originalTarget
  achievedMean=x$achievedMean
  Nindependent=x$Nindependent
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
  cat("Maximum weight was: ", max(x$weights),"\n")
  if( Nindependent >0 )
    cat("Data augmented with", Nindependent, "samples with independent variables.",
        "\nThe final weight of these samples was: ", 
        sum(weights[(length(weights)-Nindependent+1):length(weights)]), "\n")
  
  cat("----------------------------------------------------------------\n")
  
}
