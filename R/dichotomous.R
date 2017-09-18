#' @title Function dichotomous
#' @export
#' @description Simulate a dataset with
#' 5 interdependent dichotomous variables.
#' @details The output dataset has variables/columns x1 through x5.
#' All variables are dichotomous (values are either 0 or 1).
#' The marginal means are \code{probs} for x1 through x5,
#' respectively.
#' x1 and x2 are correlated such that
#' Prob(x1 = 1 and x2 = 1) is \code{prob12}.
#' x3 is nested in x1 in the sense that
#'   1. Prob(x3 = 1 | x1 = 0) = 0, and
#'   2. Prob(x3 = 1 | x1 = 1) = Prob(x3 = 1 and x1 = 1) / Prob(x1 = 1)
#'                            = Prob(x3 = 1) / Prob(x1 = 1)
#' Likewise, x4 is nested in x1 and x5 is nested in x2.
#' @param nrow number of rows in the output dataset
#' @param probs a numeric vector of 5 numbers between 0 and 1.
#' These are the marginal probabilities of x1 through x5, respectively
#' (approximate column means of the output dataset).
#' @param prob12 Prob(x1 = 1 and x2 = 1)
dichotomous <- function(
  nrow = 1e3,
  probs = c(0.5, 0.6, 0.2, 0.3, 0.4),
  prob12 = 0.4
){
  stopifnot(
    nrow > 0 &
    is.numeric(probs) &
    all(0 < probs) &
    all(1 > probs) &
    length(probs) == 5 &
    is.numeric(prob12) &
    prob12 >= 0 &
    prob12 < 1
  )

  # x1 and x2 are correlated with
  # Prob(x1 = 1 and x2 = 2) equal to prob12.
  commonprob <- matrix(c(probs[1], prob12, prob12, probs[2]), nrow = 2)
  x12 <- rmvbin(n = nrow, commonprob = commonprob)
  x1 <- x12[, 1]
  x2 <- x12[, 2]

  # x3 is nested in x1 in the sense that
  #   1. Prob(x3 = 1 | x1 = 0) = 0, and
  #   2. Prob(x3 = 1 | x1 = 1) = Prob(x3 = 1 and x1 = 1) / Prob(x1 = 1)
  #                            = Prob(x3 = 1) / Prob(x1 = 1)
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
