N <- 10000 #number of individuals
 
data.path <-  '/scratch/temporary/dtrejoba/AsynchronousGibbs'
#we generate the covariate matrix fo r the first RAR matrix
require(mvnfast)

corrMatrices <- list.files(path=data.path,full.names=T,pattern="dataR.*.RData")
print("reading list of correlation matrices")
print(corrMatrices)

create.x.matrix <- function( corr.matrix ){
	mvnfast::rmvn(N, mu = rep(0, ncol(corr.matrix)),sigma = corr.matrix , ncores=10)
}

create.b.matrix <- function( x.matrix, name=NULL){
       if(name!=NULL){
            
       }#special case for the block matrix
       else{
           
       }

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
        e <- rnorm(N,sqrt())
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
