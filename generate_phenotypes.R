N <- 1000 #number of individuals
M <- c(60,30,10 )  # effect sizes mixtures elements, 100 causals
P <- c(0.3,0.1,0.1)  # effect siyes mixtures Variance explained totaling 0.5
data.path <-  '/scratch/temporary/dtrejoba/AsynchronousGibbs/data/1k5k'
#we generate the covariate matrix fo r the first RAR matrix
require(mvnfast)

corrMatrices <- list.files(path=data.path,full.names=T,pattern="RAR|Rrandom",include.dirs=T)
corrMatrices <- corrMatrices[!grepl("_",corrMatrices)]
print("reading list of correlation matrices")
print(corrMatrices)

create.x.matrix <- function( matrix.name, corr.matrix ){
    if(grepl("random",matrix.name)){
        print("generating design matrix from the Cholesky factor")
	mvnfast::rmvn(N, mu = rep(0, ncol(corr.matrix)),sigma = corr.matrix , isChol=T, ncores=10)
    }
    else{
        print("generating design matrix from the covariance matrix")
        mvnfast::rmvn(N, mu = rep(0, ncol(corr.matrix)),sigma = corr.matrix , isChol=F, ncores=10)
    }
}

create.b.matrix <- function( x.matrix, name){
       B <- rep(0,ncol(x.matrix ))
       comps <- list()
       
       print("generating coefficients for non-block matrix")
          
       coefficients <- sample(1:ncol(x.matrix),size = sum(M) )
           
       for(i in 1: length(M)){
           tmp <- sample(coefficients,size = M[i])
           comps <- append(comps,list(tmp))
           B[tmp] <- rnorm(length(tmp),sd = sqrt(P[i]/M[i]))
           coefficients <-  setdiff(coefficients,comps[[i]])
       }
                
       list(coeff = B , comp = comps )
}

gen_lm <- function(matrixName){
	corr.matrix <- load(matrixName)
        print("generating linear model for matrix :")
        print(corr.matrix)
        print("generating design matrix X")
	X <- create.x.matrix(corr.matrix,get(corr.matrix))
        print("generating beta coefficients")
	beta <- create.b.matrix(X,corr.matrix)
        
        b <- beta$coeff
       
        print("building linear model")
        g <- scale(X) %*% b

        print("generating i.i.d noise")
        e <- rnorm(N,sd(g))

        y <- g + e
        var.y <- var(y)
        var.g <- var(g)
        var.e <- var(e)
        coeff.comp <- beta$comp
        variables <- c('X','b','e','y','var.y','var.g','var.e','coeff.comp')
        print("saving linear model")
        dest.file <- paste(data.path,paste("/",paste(corr.matrix,"_linear_model.RData",sep=""),sep = "" ),sep="")
        save(list=variables,file=dest.file)
        print("linear model saved in ")
        print(dest.file)
}


lapply(corrMatrices,gen_lm) 
