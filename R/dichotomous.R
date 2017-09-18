#' @title Function dichotomous
#' @export
#' @description Simulate a dataset with
#' 5 interdependent dichotomous variables.
#' @details The output dataset has variables/columns x1 through x5.
#' All variables are dichotomous (values are either 0 or 1).
#' The marginal means are \code{probs} for x1 through x5,
#' respectively.
#' x1 and x2 are correlated with correlation \code{cor12}.
#' x3 is nested in x1 in the sense that
#'   1. P(x3 = 1 | x1 = 0) = 0, and
#'   2. P(x3 = 1 | x1 = 1) = P(x3 = 1 and x1 = 1) / P(x1 = 1)
#'                         = P(x3 = 1) / P(x1 = 1)
#' Likewise, x4 is nested in x1 and x5 is nested in x2.
#' @param nrow number of rows in the output dataset
#' @param probs a numeric vector of 5 numbers between 0 and 1.
#' These are the approximate column means of the output dataset.
#' @param cor12 correlation between x1 and x2. Must be at least 0
#' and less than 1.
dichotomous <- function(
  nrow = 1e3,
  probs = c(0.5, 0.6, 0.2, 0.3, 0.4),
  cor12 = 0.25
){
  stopifnot(
    nrow > 0 &
    is.numeric(probs) &
    all(0 < probs) &
    all(1 > probs) &
    length(probs) == 5 &
    is.numeric(cor12) &
    cor12 >= 0 &
    cor12 < 1
  )
  
  # x1 and x2 are correlated with correlation cor12.
  prob12 <- matrix(c(probs[1], cor12, cor12, probs[2]), nrow = 2)
  x12 <- rmvbin(n = nrow, commonprob = prob12)
  x1 <- x12[, 1]
  x2 <- x12[, 2]
  
  # x3 is nested in x1 in the sense that
  #   1. P(x3 = 1 | x1 = 0) = 0
  #   2. P(x3 = 1 | x1 = 1) = P(x3 = 1 and x1 = 1) / P(x1 = 1)
  #                         = P(x3 = 1) / P(x1 = 1)
  x3 <- rep(0, nrow)
  x3[x1 == 1] <- rbinom(n = sum(x1), size = 1, prob = probs[3]/probs[1])
  
  # x4 is likewise nested in x1.
  x4 <- rep(0, nrow)
  x4[x1 == 1] <- rbinom(n = sum(x1), size = 1, prob = probs[4]/probs[1])
  
  # x5 is nested in both x2.
  x5 <- rep(0, nrow)
  x5[x2 == 1] <- rbinom(n = sum(x2), size = 1, prob = probs[5]/probs[2])
  
  data.frame(
    x1 = x1,
    x2 = x2,
    x3 = x3,
    x4 = x4,
    x5 = x5
  )
}
