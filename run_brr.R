require(tidyverse)
####functions to read data from csv and run brr
args = commandArgs(trailingOnly=TRUE)
##data.path <- paste("/scratch/temporary/dtrejoba/AsynchronousGibbs/data/",args[1],sep="")
data.path <- paste("/scratch/local/monthly/dtrejoba/AsynchronousGibbs/data/",args[1],sep="")
brr.path  <- "/home/dtrejoba/repo/BRR-build/src/brr"
root.output.path <- "/scratch/temporary/dtrejoba/AsynchronousGibbs/data/"


brr.preprocess.command <- "--preprocess --preprocess-chunks 10 --thread 2 --thread-spawned 2 "
brr.sync.command <- "--analysis-type ppbayes"
brr.async.command <- "--analysis-type asyncppbayes"

if(args[2] == "long"){
    print("running long chains")
    brr.chain.command <- paste(c("--chain-length 20000 --burn-in 10000 --thin 10 ",
                                 " --thread ", args[3],
                                 " --thread-spawned ", args[4],
                                 " --decompression-tokens ", args[5],
                                 " --analysis-concurrency ", args[6],
                                 " --S 0.001,0.01"),
                               collapse = "")

}
if(args[2] == "short"){
    print("running short chains to measure convergence")
    brr.chain.command<- paste(c( "--chain-length 4000 --burn-in 1 --thin 1 ",
                                 " --thread ", args[3],
                                 " --thread-spawned ", args[4],
                                 " --decompression-tokens ", args[5],
                                 " --analysis-concurrency ", args[6],
                                 " --S 0.001,0.01"),
                               collapse = "")

}
if(args[2] =="speed"){
    print("running short chains without data writing to test for speed")
    brr.chain.command <- paste(c("--chain-length 1000 --burn-in 1000 --thin 1",
                                 " --thread ", args[3],
                                 " --thread-spawned ", args[4],
                                 " --decompression-tokens ", args[5],
                                 " --analysis-concurrency ", args[6],
                                 " --S 0.001,0.01"),
                               collapse = "")
}

phenotype.command <- function(pheno.path){paste("--pheno",pheno.path)}
data.command <- function(data.path){paste("--data-file",data.path)}
brr.output <- function(output.path){paste("--mcmc-samples",output.path)}


brr.preprocess <- function(phen,data){
    command <- paste(c(brr.path, data.command(data), phenotype.command(phen), brr.preprocess.command),collapse= " " )
    print(command)
    system(command)
}

brr.sync <- function(phen,data,output) {
    command <- paste(c(brr.path,
                       data.command(data),
                       phenotype.command(phen),
                       brr.sync.command,
                       brr.chain.command,
                       brr.output(output)),
                     collapse = " ")
   print(command)
   system( command  )
   command
}

brr.async <- function(phen,data,output) {
    command <- paste(c(brr.path,
                       data.command(data),
                       phenotype.command(phen),
                       brr.async.command,
                       brr.chain.command,
                       brr.output(output)),collapse = " ")
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

##we save the phenotype matrix with a file for each column.
to.phen.csv <- function(phen, matrix.name) {
    for(i in 1:ncol(phen)){
         file.name <- paste(c(data.path,"/",matrix.name,"y_",as.character(i),".csvphen"),collapse="")
         print(file.name)
         write.table(x = t(as.matrix(as.numeric(scale(phen[,1])))),file = file.name,quote = F, sep = ",", row.names = F, col.names = F)
         print("saved file")
         print(file.name)
         file.name
        }

}


create.csv.files <- function(x){
   m.name <- basename(x)
   matrix.name <- strsplit(m.name,split='_')[[1]][1]
   load(x)
   to.data.csv(X,matrix.name)
   print(matrix.name)
   to.phen.csv(y,matrix.name)
}

models <- list.files(data.path,pattern = "linear",full.names=T)
library(parallel)


cl <- makeCluster(5)
clusterExport(cl,c("create.csv.files","to.data.csv","to.phen.csv","data.path"))
parLapply(cl,models,create.csv.files)

csv.files <- list.files(data.path,pattern="csv",full.names = "T")
csv.files
designs <- tibble( designs = grep("X",csv.files,value =T)) %>% mutate(matrix = basename(designs)) %>%
    separate(matrix, into = c("matrix","extension"),sep="X")
phenotypes <- tibble( pheno = grep("y_",csv.files,value=T)) %>% mutate( matrix = basename(pheno)) %>%
    separate(matrix, into = c("matrix","phenid"), sep = "y") %>% separate(phenid, into = c("phenid","type"), sep = ".c")

data.sets <- phenotypes %>% inner_join(designs)

data.sets <- data.sets %>% select(designs,pheno)
#we run preprocess
apply(data.sets,MARGIN=1,FUN=function(x){brr.preprocess(x[2],x[1])})

#we run brr in parallel for preprocessing the csv files

n.chains <- 1:4

chain.commands.sync <- apply(data.sets, MARGIN=1, FUN = function(x) {
     sapply(n.chains,function(y){
         m.name <- basename(x[2])
         m.name.split <- strsplit(m.name,split = "y")[[1]]
         m.name <- paste(c(m.name.split[1],str_extract(m.name,"_[1-9]"),"s"),collapse="")
         output <- paste(
                         c(
                             root.output.path,
                             args[1],
                             "/",
                             m.name,
                             "_sync_C",
                             as.integer(y),
                             "_",
                             args[2],
                             "_p",
                             args[3]
                         ) , Collapse = ''
                     )
         print(output)
         brr.sync(x[2],x[1],output)
     })
})

#we run brr in serie so chains do not interfere between each other
chain.commands.async <- apply(data.sets, MARGIN=1, FUN = function(x) {
     lapply(n.chains,function(y){
         m.name <- basename(x[2])
         m.name.split <- strsplit(m.name,split = "y")[[1]]
         m.name <- paste(c(m.name.split[1],str_extract(m.name,"_[1-9]"),"s"),collapse="")
         output <- paste(
             c(
                 root.output.path,
                 args[1],
                 "/",
                 m.name,
                 "_async_C",
                 as.integer(y),
                 "_",
                 args[2],
                 "_p",
                 args[3]
             ) , collapse = ''
         )
         print(output)
         brr.async(x[2],x[1],output)
     })
})
