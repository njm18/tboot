#' @method print tweights_bmr   
#' @export
print.tweights_bmr <- function(x, ...) {
  weights=x$tweights$weights
  target=x$tweights$target
  originalTarget=x$tweights$originalTarget
  achievedMean=x$tweights$achievedMean
  Nindependent=x$Nindependent
  if(is.null(originalTarget)) {
    toprint= t(cbind(achievedMean, target))
    rownames(toprint) =c("Achieved Mean", "Target Mean Posterior")
    # colnames(toprint)=colnames(dataset)
  } else {
    toprint= t(cbind(achievedMean, target, originalTarget))
    rownames(toprint) =c("Achieved Mean", "Adjusted Target Mean", "Original Target Mean Posterior")
    # colnames(toprint)=colnames(dataset)
  }
  
  
  cat("----------------------------------------------------------------\n")
  cat("Object is a 'tweights_bmr' object.\n")
  cat("Optimization was successful. The weights have a sampleing\ndistribution with means close to the attemted target:\n")
  print(toprint)
  cat("Maximum weight was: ", max(weights),"\n")
  if( Nindependent >0 )
    cat("Data augmented with", Nindependent, "samples with independent variables.",
        "\nThe final weight of these samples was: ", 
        sum(weights[(length(weights)-Nindependent+1):length(weights)]), "\n")
  
  cat("----------------------------------------------------------------\n")
  
}
