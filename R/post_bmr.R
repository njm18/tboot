#' @title Function post_bmr
#' @description Simulates the joint posterior based upon a dataset and specified marginal posterior distribution of the mean of selected variables.
#' @seealso \code{\link{tweights}}
#' @export
#' @param nsims The number of posterior simulations to draw.
#' @param tweights_bmr An object of class 'tweights_bmr' created using the 'tweights_bmr' function.

post_bmr=function(nsims, tweights_bmr) {
  if(!("tweights_bmr" %in% class(tweights_bmr)))
    stop("'tweights_bmr' must be an object of class 'tweights_bmr' from the 'tweights_bmr' function.")

  p=nrow(tweights_bmr$Csqrt)
  
  z = tweights_bmr$Csqrt %*% matrix(rnorm(p*nsims),nrow=p,ncol=nsims)
  u = apply(z, 1, function(z) list(pnorm(z))) #list to stop simplification
  sims= mapply(function(u,marginal) as.numeric(quantile(marginal, probs = u[[1]])),
             u, tweights_bmr$marginal)
  colnames(sims)=names(tweights_bmr$marginal)
  return(sims)
}
