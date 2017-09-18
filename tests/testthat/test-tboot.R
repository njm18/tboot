context("tboot")

test_that("can bootstrap with correct weights", {
  set.seed(2017)
  probs <- c(.55, .65, .25, .37, .48)
  names(probs) <- paste0("x", seq_along(probs))
  prob12 <- 0.42
  x <- dichotomous(nrow = 1e3, probs = probs, prob12 = prob12)

  target <- c(.5, .7, .3, .3, .4)
  weights <- tweights(dataset = x, target = target)
  boot <- tboot(dataset = x, weights = weights, nrow = 1e6)
  rates <- colMeans(boot)
  expect_equal(unname(rates), target, tol = 1e-2)
})
