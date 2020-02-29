context("tboot")

test_that("bootstrap weights yield correct mean", {
  set.seed(2017)
  
  #Use Iris data as example
  target=c(Sepal.Length=5.5, Sepal.Width=2.9, Petal.Length=3.4)
  
  w1=tweights(dataset = iris, target = target, distance = "klqp")
  w2=tweights(dataset = iris, target = target, distance = "euchlidean")
  w3=suppressWarnings(tweights(dataset = iris, target = target, distance = "klpq")) #calls warning with this distance
  
  boot1 <- tboot(nrow = 1e6, weights = w1)
  boot2 <- tboot(nrow = 1e6, weights = w2)
  boot3 <- tboot(nrow = 1e6, weights = w3)  
  
  expect_equal(colMeans(boot1[, names(target)]), target, tol = 1e-2)
  expect_equal(colMeans(boot2[, names(target)]), target, tol = 1e-2)
  expect_equal(colMeans(boot3[, names(target)]), target, tol = 1e-2)
  
  #Try with Nindependent
  wiris1=tweights(dataset = iris, target = target, distance = "klqp", Nindependent = 10)
  wiris2=tweights(dataset = iris, target = target, distance = "euchlidean", Nindependent = 10)
  
  boot1 <- tboot(nrow = 1e6, weights = wiris1)
  boot2 <- tboot(nrow = 1e6, weights = wiris2)
  
  expect_equal(colMeans(boot1[, names(target)]), target, tol = 1e-2)
  expect_equal(colMeans(boot2[, names(target)]), target, tol = 1e-2)
  
  
  #Use new simulated data
  dataset <- data.frame(x = rnorm(100), y = rnorm(100), z=rnorm(100))
  target=c(x=0, y=0, z=0)
  w1=tweights(dataset = dataset, target = target, distance = "klqp")
  w2=tweights(dataset = dataset, target = target, distance = "euchlidean")

  boot1 <- tboot(nrow = 1e6, weights = w1)
  boot2 <- tboot(nrow = 1e6, weights = w2)

  expect_equal(colMeans(boot1), target, tol = 1e-2)
  expect_equal(colMeans(boot2), target, tol = 1e-2)
})
