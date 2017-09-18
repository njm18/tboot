context("dichotomous")

test_that("simulated dichotomous dataset has correct properties", {
  set.seed(2017)
  tol <- 1e-2
  probs <- c(.55, .65, .25, .37, .48)
  names(probs) <- paste0("x", seq_along(probs))
  prob12 <- 0.42
  x <- dichotomous(nrow = 1e5, probs = probs, prob12 = prob12)

  # Marginal probabilities
  expect_equal(colMeans(x), probs, tol = tol)

  # Correlation of x1 and x2
  expect_equal(mean(x$x1 & x$x2), prob12, tol = tol)

  # Nested variables
  probs <- unname(probs)
  expect_equal(mean(x$x3[x$x1 == 1]), probs[3] / probs[1], tol = tol)
  expect_equal(sum(x$x3[x$x1 == 0]), 0)
  expect_equal(mean(x$x4[x$x1 == 1]), probs[4] / probs[1], tol = tol)
  expect_equal(sum(x$x4[x$x1 == 0]), 0)
  expect_equal(mean(x$x5[x$x2 == 1]), probs[5] / probs[2], tol = tol)
  expect_equal(sum(x$x5[x$x2 == 0]), 0)
})
