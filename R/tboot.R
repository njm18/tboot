#' @title Function tboot
#' @description Bootstrap \code{n} rows of \code{dataset} using
#' weights \code{p}.
#' @seealso \code{\link{tweights}}
#' @export
#' @param dataset Data frame or matrix to bootstrap.
#' @param p numeric vector of row weights for resampling.
#' Must have length \code{nrow(dataset)}
#' @param n number of rows in the new bootstrapped dataset.
tboot <- function(
  dataset,
  p = rep(1 / nrow(dataset), nrow(dataset)),
  n = nrow(dataset)) {
  if (length(p) != nrow(dataset)){
    stop("length of p must equal nrow(dataset).")
  }
  index <- sample.int(
    n = length(p),
    size = n,
    prob = p,
    replace = TRUE
  )
  dataset[index, ]
}
