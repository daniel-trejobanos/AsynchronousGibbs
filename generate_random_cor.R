
M <- 50000 #covariates
N <- 10000 #individuals

require(corrplot)
source('./generateCorrelation.R')
Rcpp::sourceCpp('generateRandomCor.cpp')
file.path <-   '/scratch/cluster/monthly/dtrejoba/AsynchronousGibbs/data/'
this.seed <- as.numeric(Sys.time())
set.seed(this.seed)
save(list = c('M','N','this.seed'), file = paste(file.path,'Rmatrices.params.RData',sep=''))

#random matrix with uniformly distributed correlations
Rrandom10k50ka1L <-  generateRandomCorL(M)
save( list = 'Rrandom10k50ka1L', file = paste(file.path,'Rrandom10k50ka1L.RData', sep = '') )

library(mvnfast)
X.Rrandom10k50ka1L <- rmvn(N,mu = rep(0,M),sigma =Rrandom10k50ka1L ,isChol = T,ncores = 20)
save( list = 'X.Rrandom10k50ka1L', file = paste(file.path,'X.Rrandom10k50ka1L.RData', sep = '') )
