N <- 10000 #number of individuals
M <- c(600,300,100 )
data.path <-  '/scratch/temporary/dtrejoba/AsynchronousGibbs'
#we generate the covariate matrix fo r the first RAR matrix
require(mvnfast)

corrMatrices <- list.files(path=data.path,full.names=T,pattern="*10k50k*")
print("reading list of correlation matrices")
print(corrMatrices)

create.x.matrix <- function( corr.matrix ){
	mvnfast::rmvn(N, mu = rep(0, ncol(corr.matrix)),sigma = corr.matrix , ncores=10)
}

create.b.matrix <- function( x.matrix, name=NULL){
       if( length(grep( "RARb" ,name) ) ==0 ){
            B <- rep(0,ncol(x.matrix ))
            coefficients <- sample(1:ncol(x.matrix),size = sum(M) )
            comps <- list()
            for(i in 1: length(M)){
                tmp <- sample(coefficients,size = M[i])
                comps <- append(comps,list(tmp))
                coefficients <-  setdiff(coefficients,comps[[i]])
            }
            
                
       }#special case for the block matrix
       else{
           load()
       }
       list(coeff = , comp = )
}

gen_lm <- function(matrixName){
	corr.matrix <- load(matrixName)
        print("generating linear model for matrix :")
        print(corr.matrix)
        print("generating design matrix X")
	X <- create.x.matrix(get(corr.matrix))
        print("generating beta coefficients")
	beta <- create.b.matrix(X,corr.matrix)
        b <- beta$coeff
        print("generating i.i.d noise")
        e <- rnorm(N,sqrt(0.5))
        print("building linear model")
        g <- X %*% b$coeff
        y <- g + e
        var.y <- var(y)
        var.g <- var(g)
        var.e <- var(e)
        coeff.comp <- b$comp
        variables <- c('X','b','e','y','var.y','var.g','var.e','coeff.comp')
        print("saving linear model")
        dest.file <- paste(data.path,paste(matrixName,"_linear_model.RData",sep=""),sep = "" )
        save(list=variables,file=dest.file)
        print("linear model saved in ")
        print(dest.file)
}


lapply(corrMatrices,gen_lm) 
