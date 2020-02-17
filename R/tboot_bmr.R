#' @title Function tboot
#' @description Bootstrap \code{nrow} rows of \code{dataset} using
#' the given row-level weights.
#' @seealso \code{\link{tweights}}
#' @export
#' @param dataset Data frame or matrix to bootstrap. Rows of the dataset must be in the 
#' same order as was used for the 'tweights' call. However the dataset may include
#' additional columns not included in the 'tweights' calll.
#' @param weights an object of class 'tweights' output from the 'tweights' function.
#' @param nrow number of rows in the new bootstrapped dataset.
#' @param fillMissingAug fill in missing augmentation with primary weights resampling.
tboot_bmr=function(nrow, bmr, tol_rel_sd=.1) {
  p=nrow(bmr$Csqrt)
  z = bmr$Csqrt %*% rnorm(p)
  u = pnorm(z)
  mu= mapply(function(u,marginal) as.numeric(quantile(marginal, probs = u)),
             u, bmr$marginal)
  names(mu)=names(bmr$marginal)

  weights=suppressWarnings(tweights(bmr$tweights$X, target=mu,
                                    distance=bmr$distance,
                                    maxit = bmr$maxit,
                                    tol=bmr$tol,
                                    warningcut=bmr$warningcut,
                                    silent=TRUE,
                                    Nindependent=bmr$Nindependent,
                                    augmentWeights=bmr$tweights$augmentWeights))
  
  if(max(abs(weights$achievedMean - mu)/bmr$marginal_sd) > tol_rel_sd){
    stop("Unable to simulate accurately.")
  }

  nweights <- length(weights$weights)
  if(nweights != (bmr$Nindependent+nrow(bmr$tweights$dataset)))
    stop("length of weights must be nrow(dataset)+Nindependent.")
  index <- sample.int(n = nweights, size = nrow, prob = weights$weights,
                      replace = TRUE)
  ret=bmr$tweights$dataset[index, ]
  tmp=mu
  names(tmp)=names(bmr$marginal)
  attr(ret, "post_bmr")=tmp
  return(ret)
}
