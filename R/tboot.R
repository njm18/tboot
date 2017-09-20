#' @title Function tboot
#' @description Bootstrap \code{nrow} rows of \code{dataset} using
#' the given row-level weights.
#' @seealso \code{\link{tweights}}
#' @export
#' @param dataset Data frame or matrix to bootstrap.
#' @param weights numeric vector of row weights for resampling.
#' Must have length \code{nrow(dataset)}
#' @param nrow number of rows in the new bootstrapped dataset.
tboot <- function(
  dataset,
  weights = rep(1 / nrow(dataset), nrow(dataset)),
  nrow = base::nrow(dataset)) {
  nweights <- length(weights)
  if (nweights != nrow(dataset)){
    stop("length of weights must be nrow(dataset).")
  }
  index <- sample.int(
    n = nweights,
    size = nrow,
    prob = weights,
    replace = TRUE
  )
  dataset[index, ]
}
