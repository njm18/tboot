#' @title Function tweights_bmr
#' @description Setup the needed prequisitesin order to prepare for bayesian marginal reconstruction (including a call to tweights). Takes as input simulations from the posterior marginal distribution of variables in a dataset.
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
#' @param augmentWeights List with weights for each variable marginal. Only has effect if Nindependent!=0. Generally this should be left as NULL. 'tweights' will fill in automatically.

tweights_bmr=function(dataset, 
                      marginal, 
                      distance="klqp",
                      maxit = 1000,
                      tol=1e-8,
                      warningcut=0.05,
                      silent=FALSE,
                      Nindependent=1,
                      augmentWeights=NULL) {
  
  if(Nindependent<1)
    stop("'tweights_bmr' requires Nindependent>0 in order to keep simulation stable.")
  
  
  marginal=lapply(marginal, sort)
  target=sapply(marginal, mean)
  marginal_sd=sapply(marginal, sd)
  w=tweights(dataset, target, distance = distance, 
             maxit = maxit, tol = tol, warningcut=0.05, silent = silent,
             Nindependent=Nindependent, augmentWeights=augmentWeights)
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
  
  augmentMeans2=sapply(names(marginal), function(nm) 
    return(crossprod( X[,nm]^2, w$augmentWeights[[nm]])))
  
  m1=w$achievedMean
  
  m2=crossprod(X[1:nrow(dataset),]*w$weights[1:nrow(dataset)],X) + 
    diag(augmentMeans2)*sum(w$weights[(nrow(X)+1):length(w$weights)])
  
  V= m2 - tcrossprod(m1)
  D=diag(1/sqrt(diag(V)))
  C=D %*% V %*% D
  
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

