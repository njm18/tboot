#' @title Function \code{tweights}
#' @description Compute edge weights \code{p}
#' such that the column means of
#' \code{tboot(dataset = dataset, p = p)}
#' equal \code{target} on average.
#' @seealso \code{\link{tboot}}
#' @export
#' @param dataset Data frame or matrix to use to find row weights.
#' @param target Numeric vector of target column means.
#' @param max_iter Maximum number of iterations for \code{optim()}.
#' @details
#' This function minimizes the Kullback-Leibler divergence
#' between two sampling schemes:
#'  (P_C) classic resampling with each observation sampled with equal prob
#'  (P_W) weighted resampling with each sample i getting a different prob (p_i)
#'       D(P_C||P_W)=sum_i log(1/n) - log(p_i)
#' The function finds p_i such that:
#'       sum_i p_i = 1
#'       dataset' p = target
#'   where dataset is a N x K matrix of clinical variables.
#'   Optimization is perfromed by utilizing lagrange multipliers
#'   Output is the optimal porbability (p)
tweights <- function(
  dataset,
  target = apply(dataset, 2, mean),
  max_iter = 1e5
  ){

  if (length(target) != ncol(dataset)){
    stop("length of target must equal ncol(dataset).")
  }

  target <- c(1, target)
  x <- as.matrix(
    cbind(
      int = 1,
      dataset
    )
  )

  heuristic <- function(lambda) {
    sum(
      (constraint_value(x, lambda) - target) ^ 2
    )
  }

  start <- rep(1, ncol(x))

  # Find the lagrange multipliers (lambda) wich satisfy the constraint closely.
  tmp <- optim(
    par = start,
    fn = heuristic
  )

  optimum_lambda <- optim(
    par = tmp$par,
    fn = heuristic,
    method = "BFGS",
    control = list(
      maxit = max_iter
    )
  )

  argmin_lagrangian(x, optimum_lambda$par)
}
