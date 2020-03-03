#' @title Function tweights_bmr
#' @description Setup the needed pre-requisites in order to prepare for bayesian marginal reconstruction (including a call to tweights). Takes as input simulations from the posterior marginal distribution of variables in a dataset.
#' @seealso \code{\link{tweights}}
#' @export
#' @param dataset Data frame or matrix to use to find row weights.
#' @param marginal Must be a named list with each element a vector of simulations of the marginal distribution of the posterior mean of data in the dataset.
#' @param distance The distance measure to minimize. Must be either 'euchlidean' or 'kl' (i.e. Kullback-Leibler). 'klqp' is recomneded.
#' @param maxit Defines the maximum number of iterations for optimizing 'kl' distance.
#' @param tol Tolerance. If the achieved mean is to far from the target (i.e. as defined by tol) an error will be thrown.
#' @param warningcut Sets the cutoff for determining when a large weight will trigger a warnint.
#' @param silent Allows silencing some messages.
#' @param Nindependent Assumes the input also includes 'Nindependent'samples with independent columns. See details.
#' @details
#' Reconstructs a correlated joint posterior from simulations from a marginal posterior.
#' Algorythm is summarized more fully in the vignettes.
#' The 'Nindependent' option augments the dataset by assuming some additional specified
#' number of patients. These pateints are assumed to made up of a random bootstrapped sample
#' from the dataset for each variable marginaly leading to indepenent variables.
#' @return 
#' An object of type \code{tweights}. This object conains the following components:
#' \describe{
#'   \item{Csqrt}{Matrix square root of the covariance.}
#'   \item{tweights}{Result from the call to tweigths.}
#'   \item{marginal}{Input marginal simulations.}
#'   \item{dataset}{Formatted dataset.}
#'   \item{target}{Attempted target.}
#'   \item{distance,maxit,tol, Nindependent, warningcut}{Inputed values to 'tweights_bmr'.}
#'   \item{Nindependent}{Inputed 'Nindependent' option.}
#'   \item{augmentWeights}{Used for 'Nindependent' option weights for each variable.}
#'   \item{weights}{tilted weights for resampling}
#'   \item{originalTarget}{Will be null if target was not changed.}
#'   \item{marginal_sd}{Standard deviation of the marginals.}
#' }
#' 

tweights_bmr=function(dataset, 
                      marginal, 
                      distance="klqp",
                      maxit = 1000,
                      tol=1e-8,
                      warningcut=0.05,
                      silent=FALSE,
                      Nindependent=1) {
  
  if(Nindependent<1)
    stop("'tweights_bmr' requires Nindependent>0 in order to keep simulation stable.")
  if(!is.list(marginal))
    stop("'marginal' must be a list.")
  if(is.null(names(marginal)))
    stop("'marginal' must be a named list.")
  
  if(is.null(colnames(dataset)))
    stop("'dataset' must have named columns starting with version 0.2.0.")
  if(any( !(names(marginal) %in% colnames(dataset))))
    stop("names of elements of 'marginal' must exist in dataset.")
  
  marginal=lapply(marginal, sort)
  target=sapply(marginal, mean)
  marginal_sd=sapply(marginal, sd)
  w=tweights(dataset, target, distance = distance, 
             maxit = maxit, tol = tol, warningcut=warningcut, silent = silent,
             Nindependent=Nindependent)
  X=w$X
  
  #double check w
  if(length(w$weights) != (Nindependent+nrow(dataset)))
    stop("length of weights must be nrow(dataset)+Nindependent.")
  if(is.null(w$augmentWeights))
    stop("Attributes of weights not set correctly for 'augmentWeights.'")
  if( !any(class(w$augmentWeights) =="list") )
    stop("'augmentWeights' must be a 'list.'")
  if( is.null(names(w$augmentWeights)) )
    stop("'augmentWeights' must be a named 'list.'")
  
  
  p_independent=sum(w$weights[(nrow(X)+1):length(w$weights)])
  normalized_notindependent_weights=w$weights[1:nrow(dataset)]/(1-p_independent)
  
  #first moment
  m1_all=w$achievedMean
  m1_independent=sapply(colnames(X), function(nm) mean(X[,nm] %*%w$augmentWeights[[nm]] ))
  
  #second moment
  m2_notindependent=crossprod(X[1:nrow(dataset),]*normalized_notindependent_weights,X)
  m2_indepentdent= tcrossprod(m1_independent)
  diag(m2_indepentdent)=sapply(colnames(X), function(nm) mean( (X[,nm]^2) %*%w$augmentWeights[[nm]] ))
  m2_all=m2_notindependent*(1-p_independent) + m2_indepentdent*p_independent
  
  #variance and correlation
  V= m2_all - tcrossprod(m1_all)
  D=diag(1/sqrt(diag(V)))
  C=D %*% V %*% D
  
  if(any(is.na(diag(D)))| any(is.infinite(diag(D)))){
    stop("Implied covariance function was two small. Most likely the assumed marginal distribution is very different from the observed data in some way.")
  }
  #Note the following code is borrowed and modified from MASS mvrnorm function
  eS <- eigen(C, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  Csqrt = eS$vectors %*% diag(sqrt(pmax(ev, 0)))
  

  ret=list(Csqrt=Csqrt, tweights=w, marginal=marginal, dataset=X,
           target=target,
           distance =distance, 
           maxit = maxit, tol =tol, 
           Nindependent=Nindependent,
           warningcut=warningcut,
           marginal_sd=marginal_sd)
  class(ret)="tweights_bmr"
  return(ret)
}

