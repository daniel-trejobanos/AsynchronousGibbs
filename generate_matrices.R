#Generating the correlation matrices
# if dimension is above 1000 then the directory created will have the name 1k
args <- commandArgs(trailingOnly=TRUE)
print(args)
print(length(args))
if(length(args) != 3)
    stop("Error: three arguments are required, number of individuals and number of covariates and number of blocks for block matrix")
M <- as.integer(args[2]) #covariates
N <- as.integer(args[1]) #individuals
B <- as.integer(args[3])  #number of blocks

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

gen.save <- function(matrix.name, ld.blocks = NULL) {
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
    if(!is.null(ld.blocks)){
        block.file <- paste(
                             c(file.path,matrix.name,"_blocks.RData"),
                             collapse = ''
                      )
        save(list = "ld.blocks", file = block.file) 
        print(block.file)
    }
}

assign.matrix <-  function(matrix.name,matrix){
    assign(matrix.name, matrix, envir = .GlobalEnv)

}


#random matrix with high correlations
tempjob1 <- function(){source('./generateCorrelation.R')
    Rcpp::sourceCpp('generateRandomCor.cpp')
    matrix.name <- paste( 'Rrandom', size.prefix, sep ='')
    assign.matrix(matrix.name, generateRandomCorL(M) )
    gen.save(matrix.name)      
}

#AR matrix with 0.5 phi 
tempjob2 <- function(){source('./generateCorrelation.R')
    matrix.name <- paste(
                          c('RAR',size.prefix,'phi0.5'),
                          collapse = ''
                   )
   assign.matrix(matrix.name,
           generateCorrelationMatrix(M,"AR", phi = 0.5)
    )
    gen.save(matrix.name)     
}
#AR matrix with high correlation phi 0.8
tempjob3 <- function(){source('./generateCorrelation.R')
     matrix.name <- paste(
                          c('RAR',size.prefix,'phi0.8'),
                          collapse = ''
                   )
    assign.matrix(matrix.name,
           generateCorrelationMatrix(M,"AR", phi = 0.8)
    )
    gen.save(matrix.name)     
}

#AR matrix with blocks
gen.delta <- function(){
       delta <- rep(1,M)
       block.lims <- sample.int(M,B)
       delta[block.lims] <- 3
       delta
}

tempjob4 <- function(){source('./generateCorrelation.R')
    matrix.name <- paste(
                          c('RARb',size.prefix,'phi0.6'),
                          collapse = ''
                   )
    delta = gen.delta()
    assign.matrix(matrix.name,
           generateCorrelationMatrix(M,"AR", phi = 0.6, delta = delta)
    )
    gen.save(matrix.name,delta)    
}

#AR matrix with blocks
tempjob5 <- function(){source('./generateCorrelation.R')
    matrix.name <- paste(
                          c('RARb',size.prefix,'phi0.8'),
                          collapse = ''
                   )
    delta = gen.delta()
    assign.matrix(matrix.name,
           generateCorrelationMatrix(M,"AR", phi = 0.8, delta = delta)
    )
    gen.save(matrix.name, delta)    
}


functionlist<- c("tempjob5","tempjob4","tempjob3","tempjob2","tempjob1")
#clusterExport(cl,c(functionlist,'M','B'),env=environment())
lapply(functionlist,function(x){eval(call(x))})

 
