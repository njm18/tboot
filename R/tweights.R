#' @title Function \code{tweights}
#' @description Returns a vector \code{p} of resampling probabilities 
#' such that the column means of \code{tboot(dataset = dataset, p = p)}
#' equals \code{target} on average.
#' @seealso \code{\link{tboot}}
#' @export
#' @param dataset Data frame or matrix to use to find row weights.
#' @param target Numeric vector of target column means.
#' @param distance The distance to minimize. Must be either 'euchlidean' or 'kl' (i.e. Kullback-Leibler).
#' @param control Controll prameter to be passed into \code{optConstr()}.
#' @param tol Tolerance. If the achieved mean is to far from the target (i.e. as defined by tol) an error will be thrown.
#' @param silent Allows silencing some messages.
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


tweights <-function(
  dataset,
  target = apply(dataset, 2, mean),
  distance="euchlidean",
  control = list(maxit=10000, reltol=1e-16),
  tol=1e-5,
  silent=FALSE) {
  oldTarget=NULL #will be replaced if we need to fix the target
  dataset=as.matrix(dataset)
  if(!is.numeric(dataset))
    stop("All columns of 'dataset' must be numeric.")
  #Check input
  if (length(target) != ncol(dataset)){
    stop("length of target must equal ncol(dataset).")
  }
  
  if(  sum(!is.finite(target))  > 0)
    stop("'target' is not valid. 'target' must be finite non missing.")
  
  if(  sum(!is.finite(dataset))  > 0)
    stop("'dataset' is not valid. 'dataset' must be finite non missing for all values.")
  
  if( any(apply(dataset, 2, var)==0) )
    stop("At least on column of 'dataset' has no variation (i.e. all subjects are the same for at least one column). Consider removing a column and its 'target'.")
  
  
  #Include probability constraint to sum to 1
  b <- c(1, target)
  A <- (as.matrix(
    cbind(
      int = 1,
      dataset
    )
  ))
  n=nrow(dataset)
  
  #Find best weights usingn euchlidian distance - if the target is infeasible try using ipop to get a closer.
  opt <- tryCatch(
    solve.QP(dvec =rep(1/n,n), Dmat =diag(n),
             Amat=cbind(A,diag(n)),
             bvec=c(b,rep(0,n)), meq =length(b), factorized = TRUE),
    error = function(e) NULL)
  
  if(is.null(opt)) {
    oldTarget=target
    target = .how_close(dataset, target)  #find a target that is achievable
    b <- c(1, target)
    opt=solve.QP(dvec =rep(1/n,n), Dmat =diag(n),
                 Amat=cbind(A,diag(n)),
                 bvec=c(b,rep(0,n)), meq =length(b), factorized = TRUE) #This time error out if you cant find it
  }
  
  
  # if(how(opt)!="converged")
  #   stop("'ipop' did not converge.")
  
  if(distance=="euchlidean") {
    return(.print_ret(weights=opt$solution, dataset, target, oldTarget, silent))
  } else if(distance=="kl") {
    
    for(i in 10:1) {
      pi=((i-1)*n*opt$solution+1)/(i*n) # slightly regularize to form a valid starting value
      lambdaStart= as.vector( solve(crossprod(A), crossprod(A , (1/pi ) ) ) )
      pi=as.vector(1/( A %*% lambdaStart))
      if(all(pi>=0))
        break
    }
    
    
    #Parameters for transformed problem to (hopefuly) make the optimization numerically stable
    s=svd(A)
    if( (s$d[length(s$d)]/s$d[1]) < 1e-6)
      warning("Matrix is ill conditioned (e.g. columns may be colinear). Please consider removing some columns.")
    vdinv = s$v %*% diag(1/s$d)
    x_star = (A %*% vdinv )
    target_star =   as.vector(t(vdinv) %*% b)
    start_star =  as.vector(solve(vdinv, lambdaStart))
    
    #Optimize constraining to the feasible space where the probability is posative
    dist=function(lambda) {
      pi=as.vector(1/(x_star %*% lambda))
      return(sqrt( sum( ( pi %*% x_star -target_star)^2) ))
    }
    optConstr = constrOptim(start_star, dist, grad=NULL, ui=x_star, ci=rep(-1,n), control=control,  outer.eps = 1e-11)
    
    
    
    if(optConstr$convergence>0 || optConstr$value >tol) {
      grad=function(lambda) {
        pi=as.vector(1/(x_star %*% lambda))
        dif=pi %*% x_star -target_star
        d= t(x_star) %*% (x_star*pi^2) 
        return(-(dif) %*% d / sqrt(sum(dif^2)))
      }
      optConstr2 = constrOptim(optConstr$par, dist, grad=grad, ui=x_star, ci=rep(-1,n), control=control, method="BFGS")
      if(optConstr2$convergence>0 || optConstr2$value >tol) {
        optConstr3 = optim(optConstr2$par, dist, gr=grad, control=control, method="BFGS")
        
        if(optConstr3$convergence>0)
          stop("Convergence not reached using 'constrOptim' or 'optim' functions.")
        if( optConstr3$value >tol)
          stop(paste("Unable to find a close enough optima. Adjust tol to avoid this error. Best tol achieved is:",  optConstr3$value, ".\n"))
      }
    }
    ret=as.vector(1/(x_star %*% optConstr$par))
    if(any(ret<0))
      stop("Optimization failed.")
    return(.print_ret(weights=ret, dataset, target, silent))
  } else
    stop("distance must be 'kl' or 'euchlidean.'")
  
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

.print_ret = function(weights, dataset, target, oldTarget, silent) {
  weights=ifelse(weights>0,weights,0)
  weights=weights/sum(weights)
  if(is.null(oldTarget)) {
    toprint= t(cbind(t(dataset) %*% weights, target))
    rownames(toprint) =c("Achieved Mean", "Target Mean")
    colnames(toprint)=colnames(dataset)
  } else {
    toprint= t(cbind(t(dataset) %*% weights, target, oldTarget))
    rownames(toprint) =c("Achieved Mean", "Adjusted Target Mean", "Original Target Mean")
    colnames(toprint)=colnames(dataset)
  }

  if(!silent){
    cat("---------------------------------------------------------\n")
    cat("Optimization was successful. The weights have a sampleing\ndistribution with means close to the attemted target:\n")
    print(toprint)
    cat("---------------------------------------------------------\n")
  }
  return(weights)
}
