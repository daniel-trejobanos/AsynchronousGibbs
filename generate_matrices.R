#Generating the correlation matrices
#
library(parallel)
library(foreach)
cl <- makeCluster(5);

M <- 50000 #covariates
N <- 10000 #individuals
B <- 4  #number of blocks
require(corrplot)
source('./generateCorrelation.R')
file.path <- '/scratch/cluster/monthly/dtrejoba/AsynchronousGibbs/data/'
this.seed <- as.numeric(Sys.time())
set.seed(this.seed)
save(list = c('M','N','this.seed'), file = paste(file.path,'Rmatrices.params.RData',sep=''))

#random matrix with high correlations
tempjob1 <- function(){source('./generateCorrelation.R',local=T)
    Rrandom10k50ka0.001 <- generateCorrelationMatrix(M,"random",alphad=0.001)
    save( list = 'Rrandom10k50ka0.001', file = paste(file.path,'Rrandom10k50ka0.001.RData', sep = '') )
    rm(Rrandom10k50ka0.001)
    gc()
    print("matrix Rrandom10k50ka0.001 generated and saved")
}
#May be to heavy to generate many matrices
#random matrix with uniformly distributed correlations
#Rrandom10k50ka1 <- generateCorrelationMatrix(M,"random",alphad=1)
#save( list = 'Rrandom10k50ka1', file = paste(file.path,'Rrandom10k50ka1.RData', sep = '') )
#rm(Rrandom10k50ka1)
#gc()

#AR matrix with 0.5 phi 
tempjob2 <- function(){source('./generateCorrelation.R',local=T)
    RAR10k50kphi0.5 <- generateCorrelationMatrix(M,"AR")
    save( list = 'RAR10k50kphi0.5', file = paste(file.path,'RAR10k50kphi0.5.RData', sep = ''))
    rm(RAR10k50kphi0.5)
    gc()
    print("matrix RAR10k50kphi0.5 generated and saved")
}
#AR matrix with high correlation
tempjob3 <- function(){source('./generateCorrelation.R',local=T)
    RAR10k50kphi0.8 <- generateCorrelationMatrix(M,"AR", phi = 0.8)
    save( list = 'RAR10k50kphi0.8', file = paste(file.path,'RAR10k50kphi0.8.RData', sep = '') )
    rm(RAR10k50kphi0.8)
    gc()
    print("matrix RAR10k50kphi0.8 generated and saved")
}

#AR matrix with blocks
tempjob4 <- function(){source('./generateCorrelation.R',local=T)
    delta <- rep(1,M)
    block.lims <- sample.int(M,B)
    delta[block.lims] <- 3
    RARb10k50kphi0.6 <- generateCorrelationMatrix(M,"AR", phi = 0.6, delta = delta)
    save( list = 'RARb10k50kphi0.6', file = paste(file.path,'RARb10k50kphi0.6.RData', sep = ''))
    rm(RARb10k50kphi0.6)
    gc()
    print("matrix RARb10k50kphi0.6 generated and saved") 
}

#AR matrix with blocks
tempjob5 <- function(){source('./generateCorrelation.R',local=T)
    delta <- rep(1,M)
    block.lims <- sample.int(M,B)
    delta[block.lims] <- 3
    RARb10k50kphi0.8 <- generateCorrelationMatrix(M,"AR", phi = 0.8, delta = delta)
    save( list = 'RARb10k50kphi0.8', file = paste(file.path,'RARb10k50kphi0.8.RData', sep = ''))
    save(list=c('block.lims','B'),file = paste(file.path,'Block_info_RARb10k50kphi0.8.RData',sep= ''))
    rm(RARb10k50kphi0.8)
    gc()
    print("matrix RARb10k50kphi0.8 generated and saved")
}

temp.list <- lapply(ls(pattern="temp"),get)

functionlist<- c("tempjob5","tempjob4","tempjob3","tempjob2","tempjob1")
clusterExport(cl,c(functionlist,'M','B'),env=environment())
lapply(functionlist,function(x){eval(call(x))})
