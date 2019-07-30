
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

#will create a matrix with nrow equals the levels of sparsity 
create.b.matrix <- function( x.matrix, name){
       B <- matrix(
           rep(
               0, ncol(x.matrix ) * length(number.b)),
           nrow = length(number.b) )
       

       comps <- list()
       
       print("generating coefficients for non-block matrix")
       
       total.number.effects <- ceiling( ncol(x.matrix) * number.b )
       effects.per.mixture <- lapply(total.number.effects, function(x) { ceiling( x * M ) } )
        ## total.number.effects contains an array, each element is the number of non-zero effects per level
        ## effects.per mixture contains a list each element is a length-3 array of non-zero elements
       ## fill up the effect matrix, by first setting the non-zero coefficients and then assigning mixture
       for(j in 1: length(number.b)){
           coefficients <- sample(1:ncol(x.matrix),size = sum(total.number.effects[j]) )
           mixture.comps <- list()
           for(i in 1: length(M)){
               tmp <- sample(coefficients,size = effects.per.mixture[[j]][i])
               mixture.comps <- append(mixture.comps,list(tmp))
               B[j,tmp] <- rnorm(length(tmp),sd = sqrt(P[i]/effects.per.mixture[[j]])) #divide the PVE by the effects in mixture
               coefficients <-  setdiff(coefficients,mixture.comps[[i]])
           }
           comps <- append(comps, mixture.comps)
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
        #we create a phenotype for each row of b
        var.y <- rep(0,nrow(beta$coeff))
        var.g <- rep(0,nrow(beta$coeff))
        var.e <- rep(0,nrow(beta$coeff))
        y <-  matrix(0, nrow = nrow(X), ncol = nrow(beta$coeff))
        b <- beta$coeff
        coeff.comp <- beta$comp
        for( i in 1: nrow(beta$coeff)){
       
            print("building linear model")
            g <- scale(X) %*% b[i,]

            print("generating i.i.d noise")
            e <- rnorm(N,sd(g))

            y[,i] <- g + e
            var.y[i] <- var(y[,i])
            var.g[i] <- var(g)
            var.e[i] <- var(e)
            
        }
       
        variables <- c('X','b','y','var.y','var.g','var.e','coeff.comp')
        print("saving linear model")
        dest.file <- paste(data.path,paste("/",paste(corr.matrix,"_linear_model.RData",sep=""),sep = "" ),sep="")
        save(list=variables,file=dest.file)
        print("linear model saved in ")
        print(dest.file)
}


number.b <- c(0.10,0.50,0.9) #effect sizes numbers as percentages of the total number of variables
M <- c(.60,.30,.10 )  # effect sizes mixtures elements, as percentages of the total number of non-zero elements 
P <- c(0.3,0.1,0.1)  # effect sizes mixtures Variance explained totaling 0.5
args <-  commandArgs(trailingOnly = TRUE)
print(args)
data.path <- paste("/scratch/temporary/dtrejoba/AsynchronousGibbs/data/",args[2],sep="")
print(length(args))
if(length(args)!=2)
    stop("Errror, wrong number of arguments: number of individuals and size of data set is need")
N <- as.integer(args[1])

    ##we generate the covariate matrix fo r the first RAR matrix
require(mvnfast)
corrMatrices <- list.files(path=data.path,full.names=T,pattern="RAR|Rrandom",include.dirs=T)
corrMatrices <- corrMatrices[!grepl("_",corrMatrices)]
corrMatrices <- grep(corrMatrices,pattern="RData",value=T)
print("reading list of correlation matrices")
print(corrMatrices)
lapply(corrMatrices,gen_lm) 

