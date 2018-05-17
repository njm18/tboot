context("edge-cases")

test_that("incorrect input to tboot", {
  dataset <- data.frame(x = rnorm(100), y = rnorm(100))
  expect_error(tboot(dataset = dataset, weights = 1))
})

test_that("incorrect input to tweights", {
  dataset <- data.frame(x = rnorm(100), y = rnorm(100))
  expect_error(tweights(dataset = mtcars, target = 1))
  expect_error(tweights(dataset = dataset, target = c(0,0), distance ="bob"))
  dataset[1,1]="NA"
  expect_error(tweights(dataset = dataset, target = c(0,0)))
})
