#Generating the correlation matrices
#
library(parallel)
library(foreach)
cl <- makeCluster(5);

M <- 5000 #covariates
N <- 1000 #individuals
B <- 50  #number of blocks

prefix.k <- function(num) {
    pre <- as.character(ceiling(num/1000))
    paste(
       pre,
       "k",
       sep = ""
    )
}

p.m <- ifelse(
          M >= 1000,
          prefix.M <- prefix.k(M),
          as.character(M)
       )  

p.n <- ifelse(
          N >= 1000,
          prefix.M <- prefix.k(N),
          as.character(N)
       )
size.prefix <- paste(p.n, p.m, sep = "")

paste(p.n,p.m,sep='')
require(corrplot)
source('./generateCorrelation.R')


file.path <- '/scratch/temporary/dtrejoba/AsynchronousGibbs/data/'
file.path <- paste(
                    c( file.path, size.prefix, "/"),
                    collapse = ""
             )
dir.create(file.path)


this.seed <- as.numeric(Sys.time())
set.seed(this.seed)
save(list = c('M','N','this.seed'), file = paste(file.path,'Rmatrices.params.RData',sep=''))

gen.save <- function(matrix.name) {
     matrix.file.name <- paste(
                               c(file.path, matrix.name, '.RData'),
                               collapse = ''
                        )
    save( list = matrix.name, file = matrix.file.name )
    print(
        paste(
            c("matrix", matrix.name, "generated and saved"),
            collapse = ' '
        )
    )
    print(matrix.file.name)
}
    

#random matrix with high correlations
tempjob1 <- function(){source('./generateCorrelation.R',local=T)
    Rcpp::sourceCpp('generateRandomCor.cpp')
    matrix.name <- paste( 'Rrandom', size.prefix, sep ='')
    assign(matrix.name, generateRandomCorL(M) )
    gen.save(matrix.name)      
}

#AR matrix with 0.5 phi 
tempjob2 <- function(){source('./generateCorrelation.R',local=T)
    matrix.name <- paste(
                          c('RAR',size.prefix,'phi0.5'),
                          collapse = ''
                   )
    assign(matrix.name,generateCorrelationMatrix(M,"AR", phi = 0.5))
    gen.save(matrix.name)     
}
#AR matrix with high correlation phi 0.8
tempjob3 <- function(){source('./generateCorrelation.R',local= T)
     matrix.name <- paste(
                          c('RAR',size.prefix,'phi0.8'),
                          collapse = ''
                   )
    assign(matrix.name,generateCorrelationMatrix(M,"AR", phi =0.8))
    gen.save(matrix.name)     
}

#AR matrix with blocks
gen.delta <- function(){
       delta <- rep(1,M)
       block.lims <- sample.int(M,B)
       delta[block.lims] <- 3
}

tempjob4 <- function(){source('./generateCorrelation.R',local=T)
   
    RARb10k50kphi0.6 <- generateCorrelationMatrix(M,"AR", phi = 0.6, delta = delta)
    save( list = 'RARb10k50kphi0.6', file = paste(file.path,'RARb10k50kphi0.6.RData', sep = ''))
    save(list=c('block.lims','B'),file = paste(file.path,'Block_info_RARb10k50kphi0.6.RData',sep= ''))    
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
#clusterExport(cl,c(functionlist,'M','B'),env=environment())
lapply(functionlist,function(x){eval(call(x))})

 #eval(call("tempjob1"))
###eval(call("tempjob4"))
#eval(call("tempjob5"))
