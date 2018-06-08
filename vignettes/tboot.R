## ----junk, cache = FALSE-------------------------------------------------
cat("Asfdasdfasfdfds")

## ----sim, cache = FALSE--------------------------------------------------
library(tboot)
set.seed(2018)
color   <- sample(c("brown", "green", "blue"), 300, replace = TRUE)
quant1  <- rnorm(300) + ifelse(color=="red", 1, 0)
quant2  <- rnorm(300) + quant1*.5
bin1    <- ifelse(quant1+rnorm(300) > 1, 1, 0)
bin2    <- ifelse(quant2+rnorm(300) > 1, 1, 0)
simData <- data.frame(color, quant1, quant2, bin1, bin2)
head(simData)

## ----dataset, cache = FALSE----------------------------------------------
dataset=as.matrix(cbind(
    colorBlue=ifelse(simData$color=="blue",1,0),
    colorBrown=ifelse(simData$color=="brown",1,0),
    simData[,-1]))
colMeans(dataset)

## ----target, cache = FALSE-----------------------------------------------
target <- rep(0.4, ncol(dataset))

## ----weights, cache = FALSE----------------------------------------------
weights <- tweights(dataset = dataset, target = target)

## ----boot, cache = FALSE-------------------------------------------------
boot <- tboot(dataset = dataset, weights = weights, nrow = 1e5)

## ----compare, cache = FALSE----------------------------------------------
colMeans(boot)

## ----hist, cache = FALSE-------------------------------------------------
hist(weights, breaks=25)
abline(v=1/300,col="red")

## ----targetweight, cache = FALSE-----------------------------------------
weights <- tweights(dataset = dataset[,c("quant2", "bin1")], 
                    target = c(0.4, 0.5))

## ----boot2, cache = FALSE------------------------------------------------
boot <- tboot(dataset = dataset, weights = weights, nrow = 1e5)
rbind("dataset mean"    = colMeans(dataset),
      "tbootstrap mean" = colMeans(boot))

