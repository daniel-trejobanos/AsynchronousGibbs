## ----eval=F--------------------------------------------------------------
## generateCorrelationMatrix(M,type,alphad=0.5,phi=0.5,sparse=F,thres=1e-4, delta = rep(1,M) )
## generateCoVarianceMatrix(corMatrix, sigmas =  rep(1,ncol(corMatrix)))

## ---- fig.cap= "correlation matrix random alphad=100"--------------------
require(corrplot)
source('./generateCorrelation.R')∑
corrplot(generateCorrelationMatrix(20,"random",alphad=100),tl.pos = "n")


## ---- fig.cap ="correlation matrix random alphad=0.001"------------------
corrplot(generateCorrelationMatrix(20,"random",alphad=0.001),tl.pos = "n")

## ----fig.cap= "Correlation matrix AR"------------------------------------
corrplot(generateCorrelationMatrix(20,"AR"),tl.pos = "n")

## ----fig.cap= "correlation matrix AR phi = 0.9"--------------------------
corrplot(generateCorrelationMatrix(20,"AR",phi = 0.9),tl.pos = "n")

## ----fig.cap= "correlation matrix AR with block structure"---------------
block_delta <- c(rep(1,4),3,rep(1,9),2,rep(1,5)) # now you have steps of size 3 at pos. 5, and 2 at pos.15
corrplot(generateCorrelationMatrix(20,"AR",phi=0.8 ,delta=block_delta),tl.pos = "n")

## ---- fig.cap= "correlation and covariance matrix with column sds =1"----
R <- generateCorrelationMatrix(20,"random",alphad=0.001)
Sigma <- generateCoVarianceMatrix(R)
corrplot(R,is.corr = F,tl.pos = "n",title = "R")
corrplot(Sigma,is.corr = F,tl.pos = "n",title = "Sigma")

## ----fig.cap=" original covariance matrix"-------------------------------
require(MASS)
M <- 20
N <- 10

R <- generateCorrelationMatrix(M,"random",alphad=0.001)
Sigma <- generateCoVarianceMatrix(R)
X <- mvrnorm(n = N, rep(0,M), Sigma)
corrplot(Sigma,is.corr = F, tl.pos = F, title = "Sigma")

## ---- fig.cap= "XtX for scaled and unscaled matrix X for N=10 and N=1000"----
par(mfrow=c(2,2))
corrplot( t(X) %*% X/N, is.corr=F, tl.pos =F, title = "un-scaled N=10")
corrplot( t(scale(X)) %*% scale(X)/N, is.corr=F, tl.pos =F, title = "scaled N=10")

N<-1000
X <- mvrnorm(n = N, rep(0,M), Sigma)
corrplot( t(X) %*% X/N, is.corr=F, tl.pos =F, title ="un-scaled N=1000")
corrplot( t(scale(X)) %*% scale(X)/N, is.corr=F, tl.pos =F, title = "scaled N=10000")

## ----fig.cap="original covariance matrix"--------------------------------
M <- 20
N <- 10
block_delta <- c(rep(1,4),3,rep(1,9),2,rep(1,5))
R <- generateCorrelationMatrix(M,"AR",phi = 0.8 ,delta = block_delta)
par(mfrow=c(1,1))
Sigma <- generateCoVarianceMatrix(R)
corrplot(Sigma,is.corr = F, tl.pos = F,title = "Sigma")


## ----fig.cap="Sigma, XtX for scaled and unscaled matrix X for N=10 and N=1000"----
X <- mvrnorm(n = N, rep(0,M), Sigma)
par(mfrow=c(2,2))

corrplot( t(X) %*% X/N, is.corr=F, tl.pos =F, title = "un-scaled N=10")
corrplot( t(scale(X)) %*% scale(X)/N, is.corr=F, tl.pos =F, title = "scaled N=10")

N<-1000
X <- mvrnorm(n = N, rep(0,M), Sigma)
corrplot( t(X) %*% X/N, is.corr=F, tl.pos =F, title ="un-scaled N=1000")
corrplot( t(scale(X)) %*% scale(X)/N, is.corr=F, tl.pos =F, title = "scaled N=10000")
dev.off()

