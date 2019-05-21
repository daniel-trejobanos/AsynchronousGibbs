
require(BayesRRcpp)#just for testing

MT = 2000
N = 5000
M = 1000
X <- matrix(rbinom(n=MT*N, 2, .5), ncol=MT)
X1<-X
X[X == 2] <- "GG"
X[X == 1] <- "AG"
X[X == 0] <- "AA"
b <- c(rnorm(M,0,sqrt(0.4/M)),rep(0,(MT-M))) # M non-zero effects
g <- scale(X1) %*% b;  
e <- rnorm(N,0, sqrt(0.6)) #residuals variance
y = g + e  #we simulate phenotype

Y<-scale(y) #we center and scale to std=1 to interpret results in scale of sd
XX <- scale(X1) #we center and scale the columns of X
sigma0=0.01# prior  variance of a zero mean gaussian prior over the mean mu NOT IMPLEMENTED YET
v0E= 7 # degrees of freedom over the inv scaled chi square prior over residuals variance
s02E = 0.6 #scale of the inv scaled chi square prior over residuals variance
v0G = 0.0001 #degrees of freedom of the inv bla bla prior over snp effects
s02G = 0.001 # scale for the samecva 
cva=as.matrix(c(0.0001,0.001,0.01)) #components variance
BayesRSamplerV2("./sim1.csv",2, 20000, 15000,5,XX, Y,sigma0,v0E,s02E,v0G,s02G,cva)
##Here I will code SBC and the metrics for the stats paper
##
BayesRSamplerV2("./sim2.csv",2, 20000, 15000,5,XX, Y,sigma0,v0E,s02E,v0G,s02G,cva)
BayesRSamplerV2("./sim3.csv",2, 20000, 15000,5,XX, Y,sigma0,v0E,s02E,v0G,s02G,cva)
BayesRSamplerV2("./sim4.csv",2, 20000, 15000,5,XX, Y,sigma0,v0E,s02E,v0G,s02G,cva)

require(coda)
require(data.table)
require(geoR)
require(MCMCpack)

C1 <- as.matrix(data.table::fread("./sim1.csv"))
C2 <- as.matrix(data.table::fread("./sim2.csv"))
C3 <- as.matrix(data.table::fread("./sim3.csv"))
C4 <- as.matrix(data.table::fread("./sim4.csv"))
#we select first the sigmas
param.family <- function( chain.matrix, family){
  if(is.character(family))
    as.mcmc(chain.matrix[,grep(family,colnames(chain.matrix))])
  else
    stop("error, family has to be character")
}
C1.sigma <- param.family(C1,"sigma")
C2.sigma <- param.family(C2,"sigma")
C3.sigma <- param.family(C3,"sigma")
C4.sigma <- param.family(C4,"sigma")

chains.sigma <- mcmc.list(C1.sigma,C2.sigma,C3.sigma,C4.sigma)

gelman.sigma <- coda::gelman.diag(chains.sigma,multivariate = T)
#This is the german plot which plots the evolution of the scale reduction factor
coda::gelman.plot(chains.sigma,autoburnin = F)

C1.beta <- param.family(C1,"beta")
C2.beta <- param.family(C2,"beta")
C3.beta <- param.family(C3,"beta")
C4.beta <- param.family(C4,"beta")


chains.beta <- mcmc.list(C1.beta,C2.beta,C3.beta,C4.beta)
#here is challenging due that we cannot do the multivariate test,  the QR decomposition fails
gelman.beta <- coda::gelman.diag(chains.beta,multivariate = F)


# This is more or less the plot we want, with all betas and their respective scale reduction factor
plot(gelman.beta$psrf[,1])

## Now lets explore autocorr functions
## 
autocor.sigma <- coda::autocorr(chains.sigma,relative = F)
autocor.beta  <- coda::autocorr.diag(chains.beta,relative = F)
#matrix with rows being the lags and columns the betas.
betas <- as.matrix(autocor.beta)

##still to do neff 
coda::effectiveSize(chains.sigma)
hist(coda::effectiveSize(chains.beta))
temp<-sbc(X,num.params = MT+2,sbc.sweeps = 10,init.thin = 1,posterior.draws = 10000,max.thin = 10,target.neff = 80000)





