####functions to read data from csv and run brr 
data.path <- "/scratch/temporary/dtrejoba/AsynchronousGibbs/data"
brr.path  <- "/home/dtrejoba/repo/BRR-build/src/brr"
output.path <- ""


brr.preprocess.command <- "--preprocess --preprocess-chunks 10 --thread 2 --thread-spawned 2 --preload"
brr.sync.command <- "--analysis-type ppbayes" 
brr.async.command <- "--analysis-type asyncppbayes"
brr.chain.command <- "--chain-length 10000 --burn-in 5000 --thin 5 --thread 12 --thread-spawned 12 --S 0.001,0.01,0.1"

phenotype.command <- function(pheno.path){paste("--pheno",pheno.path)}
data.command <- function(data.path){paste("--data-file",data.path)}
brr.output <- function(output.path){paste("--mcmc-samples",output.path)}


brr.preprocess <- function(phen,data){
    command <- paste(c(brr.path, data.command(data), phenotype.command(phen), brr.preprocess.command),collapse= " " )
    print(command)
    system(command)  
}

brr.sync <- function(phen,data,output) {
   command <- paste(c(brr.path,data.command(data),phenotype.command(phen),brr.sync.command,brr.chain.command,brr.output(output)), collapse = " ")
   print(command)
   system( command  )
   command
}
         
brr.async <- function(phen,data,output) {
   command <- paste(c(brr.path,data.command(data),phenotype.command(phen),brr.async.command,brr.chain.command,brr.output(output)),collapse = " ") 
   print(command)
   system(command)
   command
}



##functions to read linear model and save it into csv file.
to.data.csv <- function(data.matrix,matrix.name) {
    file.name <- paste(c(data.path,"/",matrix.name,"X.csv"),collapse="")
    write.table(x = t(scale(data.matrix)),file = file.name,quote = F, sep = ",", row.names = F, col.names = F)
    print("saved file")
    print(file.name)
    file.name
}

to.phen.csv <- function(phen, matrix.name) {
    file.name <- paste(c(data.path,"/",matrix.name,"y.csvphen"),collapse="")
    print(file.name)
    write.table(x = t(as.matrix(as.numeric(scale(phen)))),file = file.name,quote = F, sep = ",", row.names = F, col.names = F)
    print("saved file")
    print(file.name)
    file.name
}


models <- list.files(data.path,pattern = "linear",full.names=T)
library(parallel)


cl <- makeCluster(5)
clusterExport(cl,c("to.data.csv","to.phen.csv","data.path"))
parLapply(cl,models,function(x){
   m.name <- basename(x)
   matrix.name <- strsplit(m.name,split='_')[[1]][1]
   print(matrix.name)
   load(x)
   to.data.csv(X,matrix.name)
   to.phen.csv(y,matrix.name)
})

csv.files <- list.files(data.path,pattern="csv",full.names = "T")
csv.files
data.sets <- data.frame(designs = grep("X",csv.files,value =T), pheno = grep("y.csvphen",csv.files,value=T))
data.sets
#we run preprocess
apply(data.sets,MARGIN=1,FUN=function(x){brr.preprocess(x[2],x[1])})

#we run brr
n.chains <- 1
chain.commands.sync <- apply(data.sets, MARGIN=1, FUN = function(x) {
     sapply(n.chains,function(y){
         m.name <- basename(x[1])
         m.name <- strsplit(m.name,split = "X.csv")[[1]][1]
         output <- paste(
                         c(
                            dirname(
                               as.character(
                                 x[1]
                               )
                             )
                          , "/", m.name,"_C",as.integer(y)
                         ) , collapse = ''
                     )
         print(output)
         brr.sync(x[2],x[1],output)  
     })  
})

chain.commands.async <- apply(data.sets, MARGIN=1, FUN = function(x) {
     lapply(n.chains,function(y){
         m.name <- basename(x[1])
         m.name <- strsplit(m.name,split = "X.csv")[[1]][1]
         output <- paste(
                         c(
                            dirname(
                               as.character(
                                 x[1]
                               )
                             )
                          , "/", m.name,"_C",as.integer(y)
                         ) , collapse = ''
                     )
         print(output)
         brr.async(x[2],x[1],output)  
     })  
})


