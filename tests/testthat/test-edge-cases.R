context("edge-cases")


##############Scenario Based

test_that("incorrect input to tweights", {
  dataset <- data.frame(x = rnorm(100), y = rnorm(100))
  expect_warning(tweights(dataset = iris, target = c(Sepal.Length=5, Sepal.Width=3.6)))
  expect_error(tweights(dataset = dataset, target = c(1,0)))
  expect_error(tweights(dataset = dataset, target = c(x=0,y=0), distance ="bob"))
  dataset[1,1]=NA
  expect_error(tweights(dataset = dataset, target = c(x=0,y=0)))
  dataset[1,1]="a"
  expect_error(tweights(dataset = dataset, target = c(x=0,y=0)))
})



test_that("incorrect input to tboot", {
  dataset <- data.frame(x = rnorm(100), y = rnorm(100))
  wiris=tweights(dataset = iris, target = c(Sepal.Length=5.5, Sepal.Width=2.9))
  expect_error(tboot(nrow=100, weights = wiris, dataset = dataset))
  expect_error(tboot(nrow=dataset, weights = wiris))
  expect_error(tboot(weights = wiris))
  expect_error(tboot( wiris))
  expect_error(tboot( nrow=100, weights = dataset))
})

##############Bayesian Based

test_that("incorrect input to tweights_bmr", {
  dataset <- data.frame(x = rnorm(100), y = rnorm(100))
  dataset_marginal<-list(x=rnorm(1000, sd=.2), y=rnorm(1000,sd=.2))
  iris_marginal<-list(Sepal.Length=rnorm(1000,mean=5, sd=.2),
                      Sepal.Width=rnorm(1000,mean=3,sd=.2))
  expect_warning(tweights_bmr(dataset = iris, marginal = iris_marginal))
  expect_error(tweights_bmr(dataset = dataset, marginal = c(1,0)))
  expect_error(tweights_bmr(dataset = dataset, marginal = list(rep(1,10),rep(1,10))))
  expect_error(tweights_bmr(dataset = dataset, marginal =dataset_marginal, distance ="bob"))
  dataset[1,1]=NA
  expect_error(tweights_bmr(dataset = dataset, marginal =dataset_marginal))
  dataset[1,1]="a"
  expect_error(tweights_bmr(dataset = dataset, marginal =dataset_marginal))
})

test_that("incorrect input to post_bmr", {
  iris_marginal<-list(Sepal.Length=rnorm(1000,mean=5.8, sd=.2),
                      Sepal.Width=rnorm(1000,mean=3,sd=.2))
  wiris=tweights_bmr(dataset = iris,
                     marginal = iris_marginal)
  
  expect_error(post_bmr(nsims = "asdfa", weights_bmr = wiris))
  expect_error(post_bmr(nsims = c(3,5), weights_bmr = wiris))
  expect_error(post_bmr(nsims = 10, weights_bmr = iris_marginal))
})

test_that("incorrect input to tboot_bmr", {
  dataset <- data.frame(x = rnorm(100), y = rnorm(100))
  iris_marginal=list(Sepal.Length=rnorm(1000,mean=5.8, sd=.2), Sepal.Width=rnorm(1000,mean=3,sd=.2))
  wiris=tweights_bmr(dataset = iris,
                     marginal = iris_marginal)
  expect_error(tboot_bmr(nrow="100", weights_bmr = wiris))
  expect_error(tboot_bmr(nrow=c(1,10), weights_bmr = wiris))
  expect_error(tboot_bmr(nrow=100, weights_bmr = dataset))
  expect_error(tboot_bmr(nrow=dataset, weights = wiris))
  expect_error(tboot_bmr(weights_bmr = wiris))
  expect_error(tboot_bmr( wiris))
})






