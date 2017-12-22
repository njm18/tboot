context("tboot")

test_that("can bootstrap with correct weights", {
  set.seed(2017)
  probs <- c(.55, .65, .25, .37, .48)
  names(probs) <- paste0("x", seq_along(probs))
  prob12 <- 0.42
  x <- dichotomous(nrow = 1e3, probs = probs, prob12 = prob12)

  target <- c(.5, .7, .3, .3, .4)
  weights_eu <- tweights(dataset = x, target = target)
  weights_kl <- tweights(dataset = x, target = target)
  
  boot_eu <- tboot(dataset = x, weights = weights_eu, nrow = 1e6)
  boot_kl <- tboot(dataset = x, weights = weights_kl, nrow = 1e6)
  
  rates_eu <- colMeans(boot_eu)
  rates_kl <- colMeans(boot_eu)
  
  expect_equal(unname(rates_eu), target, tol = 1e-2)
  expect_equal(unname(rates_kl), target, tol = 1e-2)
})
