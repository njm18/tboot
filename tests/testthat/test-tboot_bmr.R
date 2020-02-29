context("tboot")

test_that("bmr bootstrap weights yield correct distribution", {
  set.seed(2020)
  
  #Use Iris data as example
  marginal<-list(Sepal.Length=rnorm(10000,mean=mean(iris$Sepal.Length), sd=.4),
                 Sepal.Width=rnorm(10000,mean=mean(iris$Sepal.Width),sd=.4),
                 Petal.Length=rnorm(10000,mean=mean(iris$Petal.Length),sd=.4)
  )
  
  #Check thaat Csqrt is getting calculated correctly
  #note that the weights should not be tilted give that the mean of marginal is the same as the data
  #thus the correlations from data and from tweights should be close but not the same due to Nindependent
  #option
  w1=tweights_bmr(dataset = iris, marginal = marginal,
                  distance = "klqp",Nindependent = 10)
  calculateCorrTweights = as.vector(tcrossprod(w1$Csqrt))
  corrData=as.vector(cor(iris[,names(marginal)]))
  corrBoot=as.vector(cor(tboot(1e6,w1$tweights)[,names(marginal)]))
  expect_equal(calculateCorrTweights, corrData, tol = .05)
  expect_equal(calculateCorrTweights, corrBoot, tol = 5e-3)
  
  #Check thaat Csqrt again with different Nindependent option
  w1=tweights_bmr(dataset = iris, marginal = marginal,
                  distance = "klqp",Nindependent = 1)
  calculateCorrTweights = as.vector(tcrossprod(w1$Csqrt))
  corrData=as.vector(cor(iris[,names(marginal)]))
  corrBoot=as.vector(cor(tboot(1e6,w1$tweights)[,names(marginal)]))
  expect_equal(calculateCorrTweights, corrData, tol = .05)
  expect_equal(calculateCorrTweights, corrBoot, tol = 5e-3)
  

  #Use winsorized marginal to keep marginal within feasible region
  winsor=function(marginalSims,y)  {
    l=min(y)
    u=max(y)
    ifelse(marginalSims<l,l,ifelse(marginalSims>u,u, marginalSims))
  }
  
  marginal<-list(Sepal.Length=winsor(rnorm(10000,mean=5.8, sd=.2),iris$Sepal.Length),
                 Sepal.Width=winsor(rnorm(10000,mean=3,sd=.2), iris$Sepal.Width),
                 Petal.Length=winsor(rnorm(10000,mean=3.7,sd=.2), iris$Petal.Length)
  )
  w1=tweights_bmr(dataset = iris, marginal = marginal, distance = "klqp")
  w2=tweights_bmr(dataset = iris, marginal = marginal, distance = "euchlidean")

  post1 <- post_bmr(1e6, weights = w1)
  post2 <- post_bmr(1e6, weights = w2)
  
  #check first moment
  margin_mean=sapply(marginal, function(x) mean(x))
  expect_equal(colMeans(post1[, names(marginal)]), margin_mean, tol = 5e-3)
  expect_equal(colMeans(post2[, names(marginal)]), margin_mean, tol = 5e-3)
  
  #checking additional centered marginal moments of the marginal distribution
  momentroot=function(mat, m) apply(mat, 2, function(x) abs(mean((x-mean(x))^m))^(1/m))
  for(m in 2:4) { 
    margin_rootmoment=sapply(marginal, function(x) abs(mean((x-mean(x))^m))^(1/m))
    expect_equal(momentroot(post1[, names(marginal)],m), margin_rootmoment, tol = 5e-3)
    expect_equal(momentroot(post2[, names(marginal)],m), margin_rootmoment, tol = 5e-3)
  }
  
  #Check that the correlation of the posterior is approximately as expected
  calculateCorrTweights1 = as.vector(tcrossprod(w1$Csqrt))
  corrPost1=as.vector(cor(post1))
  calculateCorrTweights2 = as.vector(tcrossprod(w2$Csqrt))
  corrPost2=as.vector(cor(post2))
  expect_equal(calculateCorrTweights1, corrPost1, tol = 5e-3)
  expect_equal(calculateCorrTweights2, corrPost2, tol = 5e-3)

  #Check that the attributes from tboot_bmr are approximately as expected
  myattr=replicate(2e3, 
          {
            b=tboot_bmr(2,w1, tol_rel_sd = 3)
            attr(b, "post_bmr")
          })

  expect_equal(rowMeans(myattr), margin_mean, tol = 5e-3)
  m=2
  margin_rootmoment=sapply(marginal, function(x) abs(mean((x-mean(x))^m))^(1/m))
  expect_equal(momentroot(t(myattr),m), margin_rootmoment, tol = 5e-3) #didn't sim enought to get high tol on higer moments
  
  #Check that the the bootstrap goes to attribute
  boot1=tboot_bmr(1e6,w1)
  boot2=tboot_bmr(1e6,w2)
  expect_equal(attr(boot1, "post_bmr"), colMeans(boot1), tol = 5e-3)
  expect_equal(attr(boot2, "post_bmr"), colMeans(boot2), tol = 5e-3)

})
