context("edge-cases")

test_that("incorrect input to tboot", {
  dataset <- data.frame(x = rnorm(100), y = rnorm(100))
  expect_error(tboot(dataset = dataset, weights = 1))
})

test_that("incorrect input to tweights", {
  expect_error(tweights(dataset = mtcars, target = 1))
})
