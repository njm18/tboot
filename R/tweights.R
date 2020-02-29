#' @title Function \code{tweights}
#' @description Returns a vector \code{p} of resampling probabilities 
#' such that the column means of \code{tboot(dataset = dataset, p = p)}
#' equals \code{target} on average.
#' @seealso \code{\link{tboot}}
#' @export
#' @param dataset Data frame or matrix to use to find row weights.
#' @param target Numeric vector of target column means. If the 'target' is named, then all elements of names(target) should be in the dataset.
#' @param distance The distance to minimize. Must be either 'euchlidean,' 'klqp' or 'klpq' (i.e. Kullback-Leibler). 'klqp' which is expontential tilting is recomneded.
#' @param maxit Defines the maximum number of iterations for optimizing 'kl' distance.
#' @param tol Tolerance. If the achieved mean is to far from the target (i.e. as defined by tol) an error will be thrown.
#' @param warningcut Sets the cutoff for determining when a large weight will trigger a warnint.
#' @param silent Allows silencing some messages.
#' @param Nindependent Assumes the input also includes 'Nindependent'samples with independent columns. See details.
#' @details
#' Let \eqn{p_i = 1/n} be  probability of sampling subject \eqn{i} from a dataset with \eqn{n} individuals (i.e. rows of the dataset) in the classic resampling with replacement scheme.
#' Also, let \eqn{q_i} be the probability of sampling subject \eqn{i} from a dataset with \eqn{n} individuals in our new resampling scheme. Let \eqn{d(q,p)} represent a distance between the two resampling schemes.  The \code{tweights}
#' function seeks to solve the problem: 
#' \deqn{q = argmin_p d(q,p)}
#' Subject to the constraint that:
#' \deqn{ sum_i q_i = 1} and
#' \deqn{  dataset' q = target}
#' where dataset is a n x K matrix of variables input to the function.
#' 
#'   \deqn{d_euclidian(q,p) = sqrt( sum_i (p_i-q_i)^2 )}
#'   \deqn{d_kl(q,p) = sum_i (log(p_i) - log(q_i))}
#'
#' Optimization for euclidean distance is a quadratic program and utilizes the ipop function in kernLab.
#' The euclidean based solution helps form a starting value which is used along with the constOptim function 
#' and lagrange multipliers to solve the Kullback-Leibler distance optimization.
#'   Output is the optimal porbability (p)
#' 
#' The 'Nindependent' option augments the dataset by assuming some additional specified
#' number of patients. These pateints are assumed to made up of a random bootstrapped sample
#' from the dataset for each variable marginaly leading to indepenent variables. 
#' @return 
#' An object of type \code{tweights}. This object conains the following components:
#' \describe{
#'   \item{weights}{tilted weights for resampling}
#'   \item{originalTarget}{Will be null if target was not changed.}
#'   \item{target}{Actual target that was attempted.}
#'   \item{achievedMean}{Achieved mean from tilting.}
#'   \item{dataset}{Inputed dataset.}
#'   \item{X}{Reformated dataset.}
#'   \item{Nindependent}{Inputed 'Nindependent' option.}
#'   \item{augmentWeights}{Used for 'Nindependent' option weights for each variable.}
#' }
#' 
#' 
#' 
#' 
tweights <-function(
  dataset,
  target = apply(dataset, 2, mean),
  distance="klqp",
  maxit = 1000,
  tol=1e-8,
  warningcut=0.05,
  silent=FALSE,
  Nindependent=0
  ) {
  
  originalTarget=NULL #will be replaced if we need to fix the target
  originalDataset=dataset
  

  
  if(is.null(colnames(dataset)))
    stop("'dataset' must have named columns starting with version 1.")


  #Check input
  if(!is.numeric(target)) {
    stop("'target' must be a named numeric vector.")
  } else if(is.null(names(target))) {
      stop("'target' must be a named vector starting with version 1.")
  } else if(!all(names(target) %in% colnames(dataset))) {
      stop("Some names of 'target' have no match in colnames  'dataset.'")
  }
  
  dataset=dataset[,names(target)]
  dataset=as.matrix(dataset)
  
  if(!is.numeric(dataset))
    stop("All targeted columns of 'dataset' must be numeric.")

  if(  sum(!is.finite(target))  > 0)
    stop("'target' is not valid. 'target' must be finite non missing.")
  
  if(  sum(!is.finite(dataset))  > 0)
    stop("'dataset' is not valid. 'dataset' must be finite non missing for all values.")
  
  if( any(apply(dataset, 2, var)==0) )
    stop("At least on column of 'dataset' has no variation (i.e. all subjects are the same for at least one column). Consider removing a column and its 'target'.")
  
  if(Nindependent!=0) {
    if(floor(Nindependent)!=Nindependent)
      stop("'Nindependent' must be an integer.")
    # if(is.null(augmentWeights)) {
      augmentWeights=lapply(names(target), function(nm){
        ret=tweights(
          dataset[,nm, drop=FALSE],
          target = target[nm],
          distance=distance,
          maxit = maxit,
          tol=tol,
          warningcut=1,#dont warn
          silent=TRUE)#dont talk - remember it is the final check that matters
        return(ret$weights)
      })
      names(augmentWeights)=colnames(dataset)
    # } else{
    #   #augment weights was sent in so check its validity
    #   if( !any(class(augmentWeights) =="list") )
    #     stop("'augmentWeights' must be a 'list.'")
    #   if(is.null(names(target)))
    #     stop("'augmentWeights' must be a named list.")
    #   if(!all(names(target) %in% names(augmentWeights)))
    #     stop("Some names of 'target' have no match to a name in 'augmentWeights.'")
    #   augmentWeights=augmentWeights[names(target)]
    # }
    
    augmentMeans=sapply(names(target), function(nm) 
      return(crossprod( dataset[,nm], augmentWeights[[nm]])))


    
    augmentMeansrep=do.call(rbind, 
                         replicate(Nindependent, augmentMeans,
                                   simplify = FALSE))
    dataset=rbind(dataset,augmentMeansrep) 
  } else
    augmentWeights=NULL
     
  #Include probability constraint to sum to 1
  b <- c(1, target)
  A <- (as.matrix(
    cbind(
      int = 1,
      dataset
    )
  ))
  n=nrow(dataset)
  
  if(distance=="euchlidean") {
    #Find best weights using euchlidian distance - if the target is infeasible try using ipop to get a closer.
    opt <- tryCatch(
      solve.QP(dvec =rep(1/n,n), Dmat =diag(n),
               Amat=cbind(A,diag(n)),
               bvec=c(b,rep(0,n)), meq =length(b), factorized = TRUE),
      error = function(e) NULL)
    
    if(is.null(opt)) {
      originalTarget=target
      target = .how_close(dataset, target)  #find a target that is achievable
      b <- c(1, target)
      opt=solve.QP(dvec =rep(1/n,n), Dmat =diag(n),
                   Amat=cbind(A,diag(n)),
                   bvec=c(b,rep(0,n)), meq =length(b), factorized = TRUE) #This time error out if you cant find it
    }
    
    return(.print_ret(originalDataset,weights=opt$solution, dataset, target, originalTarget,
                      warningcut, silent, Nindependent, augmentWeights))
    
  } else if(distance=="klpq") {
    warning("The 'klpq' distance is difficult to optimize and may lead to unstable results. It has not been validated and is not recomended.")
    opt=.newtonKLpq(A,b,maxit,tol)
    if(opt$steps==maxit) {
      originalTarget=target
      target = .how_close(dataset, target)  #find a target that is achievable
      b <- c(1, target)
      opt=.newtonKLpq(A,b,maxit,tol)
      if(opt$steps==maxit)
        stop("Optimization failed. Maximum iterations reached.")
    }
    
    return(.print_ret(originalDataset,weights=opt$weights, dataset, target,originalTarget,
                      warningcut,silent, Nindependent, augmentWeights))
  } else if(distance=="klqp") {
    opt=.newtonKLqp(A,b,maxit,tol)
    if(opt$steps==maxit) {
      originalTarget=target
      target = .how_close(dataset, target)  #find a target that is achievable
      b <- c(1, target)
      opt=.newtonKLqp(A,b,maxit,tol)
      if(opt$steps==maxit)
        stop("Optimization failed. Maximum iterations reached.")
    }
    return(.print_ret(originalDataset,weights=opt$weights, dataset, target,originalTarget, 
                      warningcut,silent, Nindependent, augmentWeights))
  } else
    stop("distance must be 'klqp', 'klpq' or 'euchlidean.'")
}



.how_close <- function(dataset, target) {
  xstar=scale(dataset, center = FALSE)
  scl=attr(xstar,"scaled:scale")
  targetstar=target/scl
  
  #optimize and check convergence
  opt = ipop(c=-as.vector(xstar %*% targetstar), H=t(xstar),
             A=matrix(1,1,nrow(xstar)), b=1, r=0,
             l=rep(0,nrow(xstar)), u=rep(1,nrow(xstar)), maxiter=100)
  if(how(opt)!="converged")
    stop("'ipop' did not converge.")
  
  best=as.vector(t(xstar) %*% primal(opt))*scl
  warning(paste("Target apears to not be achievable. Replacing with the nearest achievable target in terms of scaled euclidian distance. New target is: \n", paste(best,collapse=", "),"\n"))
  return(best)
}

.newtonKLpq = function(A, b, maxit, tol) {
  #Parameters for transformed problem to (hopefuly) make the optimization numerically stable
  s=svd(A)
  if( (s$d[length(s$d)]/s$d[1]) < 1e-6)
    warning("Matrix is ill conditioned (e.g. columns may be colinear). Please consider removing some columns.")
  vdinv = s$v %*% diag(1/s$d)
  x_star = (A %*% vdinv )
  target_star =   as.vector(t(vdinv) %*% b)
  lambda_n = .normalize(x_star, 
                        as.vector( t(s$v %*% diag(s$d)) %*% c(1,rep(0, length(target_star)-1)) ))
  
  for(steps in 1:maxit) {
    xlambda_n=x_star %*% lambda_n      
    pi_n=as.vector(1/xlambda_n)
    tmp=x_star * (pi_n)
    negF_deriv=  crossprod(tmp)
    dif=as.vector(pi_n %*% x_star -target_star)
    if(any(pi_n<0))
      stop("Error in optimization step. Target may be to far from data.")
    
    if(any(is.na(dif)))
      stop("Error in optimization step. Target may be to far from data.")

    if(max(abs(dif))<tol)
      break
    lambda_proposal= lambda_n + .solveTrap(negF_deriv, dif)
    xlambda_proposal=x_star %*% lambda_proposal 
    boundary=xlambda_n/(xlambda_n-xlambda_proposal)
    boundary=boundary[boundary>0]
    if(length(boundary)>0) {
      boundary=min(boundary)
      if(boundary>1) {
        alpha=1
      } else {
        alpha=.5*boundary
      }
    }
    else
      alpha=1
    # if(alpha< 0 )
    #   browser()
    lambda_n=as.vector(alpha*lambda_proposal + (1-alpha)*lambda_n)
  }
  xlambda_n=x_star %*% lambda_n      
  pi_n=as.vector(1/xlambda_n)
  return(list(steps=steps, weights=pi_n))
} 



.solveTrap=function(A,b) {
  tryCatch(solve(A,b), error = function(e) {
    warning("Matrix was not invertible. Adding a small constant to the diagonal.")
    e=eigen(A, symmetric=TRUE)
    const = max(e$values)*.001
    return(e$vectors %*% (diag(1/(e$values+const)) %*% (t(e$vector) %*% b)))
  })
}


.newtonKLqp = function(A, b, maxit, tol) {
  
  #Parameters for transformed problem to (hopefuly) make the optimization numerically stable
  s=svd(A)
  if( (s$d[length(s$d)]/s$d[1]) < 1e-6)
    warning("Matrix is ill conditioned (e.g. columns may be colinear). Please consider removing some columns.")
  vdinv = s$v %*% diag(1/s$d)
  x_star = (A %*% vdinv )
  target_star =   as.vector(t(vdinv) %*% b)
  lambda_n =as.vector( t(s$v %*% diag(s$d)) %*% c(-log(nrow(A)),rep(0, length(target_star)-1)) )
  
  for(steps in 1:maxit) {
    xlambda_n=x_star %*% lambda_n      
    pi_n=as.vector(exp(xlambda_n))
    tmp=x_star * (pi_n)
    F_deriv=  crossprod(tmp, x_star)
    dif=as.vector(pi_n %*% x_star -target_star)
    if(any(is.na(dif)))
      stop("Error in optimization step. Target may be to far from data.")
    if(max(abs(dif))<tol)
      break
    lambda_n = lambda_n - .solveTrap(F_deriv, dif)
  }
  return(list(steps=steps, weights=pi_n))
} 


.print_ret = function(originalDataset,weights, dataset, target, originalTarget, 
                      warningcut, silent,
                      Nindependent, augmentWeights) {
  weights=ifelse(weights>0,weights,0) #Just in case of numeric issues
  weights=weights/sum(weights)
  achievedMean=as.vector(t(dataset) %*% weights)
  names(achievedMean)=colnames(dataset)
  if(is.null(originalTarget)) {
    toprint= t(cbind(achievedMean, target))
    rownames(toprint) =c("Achieved Mean", "Target Mean")
    colnames(toprint)=colnames(dataset)
  } else {
    toprint= t(cbind(achievedMean, target, originalTarget))
    rownames(toprint) =c("Achieved Mean", "Adjusted Target Mean", "Original Target Mean")
    colnames(toprint)=colnames(dataset)
  }

  if(!silent){
    cat("----------------------------------------------------------------\n")
    cat("Optimization was successful. The weights have a sampleing\ndistribution with means close to the attemted target:\n")
    print(toprint)
    cat("Maximum weight was: ", max(weights),"\n")
    if( Nindependent >0 )
      cat("Data augmented with", Nindependent, "sample(s) with independent variables.",
          "\nThe final weight of the indpendent sample(s) was: ", 
          sum(weights[(length(weights)-Nindependent+1):length(weights)]), "\n")
    
    cat("----------------------------------------------------------------\n")
  }
  if(any(weights>warningcut))
    warning(paste0("Some of the weights are larger than ", warningcut, 
                   ". Thus your bootstrap sample may be overly dependent on a few samples. See vignette."))

  
  ret = list(weights=weights,
             target=target,
             dataset=originalDataset,
             X=dataset[1:nrow(originalDataset),],
             originalTarget=originalTarget,
             achievedMean=achievedMean,
             Nindependent = Nindependent,
             augmentWeights = augmentWeights) 
  class(ret)="tweights"
  return(ret)
}

.normalize=function(X,lambda) {
  pi=1 / (X %*% lambda)
  const=sum(pi)
  return(lambda*const)
}

