require(clusterGeneration)

generateCorrelationMatrix <- function(M,type,alphad=0.5,phi=0.5,sparse=F,thres=1e-4, delta = rep(1,M) )
{
  R=matrix(rep(0,M*M),nrow = M,ncol = M)
  if(type=="random")
  {
    R<-clusterGeneration::rcorrmatrix(M,alphad)
  }
 if(type=="AR")
 #AR process
 {
  #delta between timepoints
  sigma <- 1
  t <- cumsum(delta)
  Sigma <- sigma^2/(1-phi^2)*phi^abs(outer(t,t,"-"))
  Dinv <- 1/sqrt(diag(Sigma))
  R <- diag(Dinv) %*% Sigma %*% diag(Dinv)
 }
  if(sparse)
  {
    R[R<thres]<-0
  }  
  as.matrix(R)
}

generateCoVarianceMatrix <- function(corMatrix, sigmas =  rep(1,ncol(corMatrix)))
{
  diag(sigmas) %*% corMatrix %*% diag(sigmas) 
}