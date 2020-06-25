#' @title Function tboot_bmr
#' @description Bootstrap \code{nrow} rows of \code{dataset} using
#' the given row-level weights.
#' @seealso \code{\link{tweights}}
#' @export
#' @param nrow number of rows in the new bootstrapped dataset.
#' @param weights_bmr an object of class 'tweights' output from the 'tweights' function.
#' @param tol_rel_sd An error will be called if for some simulation if the target is not achievable with the data. However, the error will only be called if max absolute difference releative to the marginal standard is greater than specified.
#' @details
#' Simulates a dataset by first simulating from the posterior distribution of the column means and then simulating a dataset with that underlying mean. Details a further documented in the vignette.
#' @return 
#' A simulated dataset with 'nrow' rows. The underlying 'true' posterior parameter value is an attribute which can be extracted useing \code{attr(ret, "post_bmr")} where 'ret' is the matrix.

tboot_bmr=function(nrow, weights_bmr, tol_rel_sd=.05) {
  if(length(tol_rel_sd)!=1 | !is.numeric(tol_rel_sd))
    stop("'tol_rel_sd' must be a numeric fector")
  if(missing(nrow))
    stop("'nrow' is missing")
  if(missing(weights_bmr))
    stop("'weights_bmr' is missing")
  if(!("tweights_bmr" %in% class(weights_bmr)))
    stop("'weights_bmr' must be an object of class 'tweights_bmr' from the 'tweights_bmr' function.")
  if(!is.numeric(nrow))
    stop("'nrow' must be numeric.")
  if(length(nrow)!=1)
    stop("'nrow' must be length 1.")
  
  p=nrow(weights_bmr$Csqrt)
  z = weights_bmr$Csqrt %*% rnorm(p)
  u = pnorm(z)
  mu= mapply(function(u,marginal) as.numeric(quantile(marginal, probs = u)),
             u, weights_bmr$marginal)
  names(mu)=names(weights_bmr$marginal)

  weights=suppressWarnings(tweights(weights_bmr$tweights$X, target=mu,
                                    distance=weights_bmr$distance,
                                    maxit = weights_bmr$maxit,
                                    tol=weights_bmr$tol,
                                    warningcut=weights_bmr$warningcut,
                                    silent=TRUE,
                                    Nindependent=weights_bmr$Nindependent))
  
  if(max(abs(weights$achievedMean - mu)/weights_bmr$marginal_sd) > tol_rel_sd){
    stop("Unable to simulate accurately.")
  }
  
  ret=tboot(nrow=nrow,
            weights=weights,
            dataset=weights$X, #for now don't include extra columns
            fillMissingAug=FALSE)
  tmp=mu
  names(tmp)=names(weights_bmr$marginal)
  attr(ret, "post_bmr")=tmp
  return(ret)
}
