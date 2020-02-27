#' @export
print.tweights_bmr <- function(obj) {
#browser()
  weights=obj$tweights$weights
  target=obj$tweights$target
  originalTarget=obj$tweights$originalTarget
  achievedMean=obj$tweights$achievedMean
  Nindependent=obj$Nindependent
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
