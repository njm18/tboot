#' Twisted Bootstrap.
#' Weighted resampling of data rows
#' to achieve target column means.
#' @docType package
#' @name tboot-package
#' @author William Michael Landau \email{will.landau@@lilly.com}
#' @author Nathan Morris \email{morris_nathan@lilly.com}
#' @references \url{https://github.com/wlandau-lilly/tboot}
#' @importFrom stats optim rbinom var pnorm quantile rnorm sd
#' @importFrom kernlab ipop primal how
#' @importFrom quadprog solve.QP
#' @suggests quadprog solve.QP
NULL