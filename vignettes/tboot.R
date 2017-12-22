## ----setup, cache = TRUE-------------------------------------------------
library(tboot)
set.seed(2017)

## ----sim, cache = TRUE---------------------------------------------------
dataset <- dichotomous(nrow = 200)

## ----dataset, cache = TRUE-----------------------------------------------
head(dataset)
colMeans(dataset)

## ----target, cache = TRUE------------------------------------------------
target <- rep(0.4, ncol(dataset))

## ----weights, cache = TRUE-----------------------------------------------
weights <- tweights(dataset = dataset, target = target)

## ----boot, cache = TRUE--------------------------------------------------
boot <- tboot(dataset = dataset, weights = weights, nrow = 1e5)

## ----compare, cache = TRUE-----------------------------------------------
colMeans(boot)

## ----hist, cache = TRUE--------------------------------------------------
hist(weights, breaks=25)
abline(v=1/200,col="red")

