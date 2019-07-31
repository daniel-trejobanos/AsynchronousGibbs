#data.path <- "/scratch/temporary/dtrejoba/AsynchronousGibbs/data/1k5k"
data.path <- "/scratch/temporary/dtrejoba/AsynchronousGibbs/1k5k_0107"
files.path <- list.files(data.path, pattern = "_C", full.names = T)
matrices <-   list.files(data.path, pattern = "_C")

## we build a data frame where we hold the matrix name and the parameters and the type


require(coda)
require(tidyverse)
#source('./sbc.R')

chains <- data.frame(directory = files.path, matrices = matrices) %>%
           separate( "matrices", into = c("matrix","method","chain"), sep="_") %>%
           separate( "matrix", into = c("matrix", "phi"), sep = "phi") %>%
           separate( "matrix", into = c("type","variables"), sep ="1k") %>%
           mutate_at(vars(directory),as.character)
#function to extract a family of parameters from a mcmc chain
param.family <- function( chain.matrix, family){
#  if(is.character(family))
   as.mcmc(select(chain.matrix,matches(family)))
 # else
  #  stop("error, family has to be character")
}
posterior.path <- paste(data.path,"/post",sep ="")
dir.create(posterior.path)
#we compute gelman R for each matrix
chains.gelman <- function(mcmc.chains,variable.family, multivariate = F){
                  mcmc.chains %>%
                  group_by(type,phi,method) %>%
                  do({ function(y){
                                  tmp <- coda::gelman.diag( mcmc.list(lapply(y, function(x){
                                                     param.family(data.table::fread(x),variable.family)
                                                 })),
                                                autoburnin = F,
                                                multivariate = multivariate)
                                   tmp.df <- as.data.frame(tmp$psrf)
                                   tmp.df$var <- rownames(tmp$psrf)
                                   tmp.df$"Upper C.I." <- NULL
                                   tmp.df %>% spread(var, "Point est.")

                                }
                               } (.$directory)
                  )



                 }

#we compute the effective sample size for each
tmp<-chains.gelman(chains,"sigma|mu")

chains.effective <-  function(mcmc.chains,variable.family){
                  tmp <- lapply(mcmc.chains$directory, FUN = function(x) { coda::effectiveSize(param.family(data.table::fread(x),variable.family))})
                  cbind(chains,do.call(rbind,tmp))
                 }

chains.effective(chains, "sigma|mu")

avg.effective <- function(mcmc.chains, variable.family){
                  chains.effective(mcmc.chains, variable.family) %>%
                     group_by(type,phi,method) %>%
                      summarise_if(is.numeric , c(mean,sd))
    }

tmp <-  avg.effective(chains, "sigma|mu")

# Here we want to read the chains and thin them and compute the parameters summaries
thin.and.summ <- function(mcmc.chains,variable.family,length.out=1001){
         mcmc.chains %>%
             group_by(type,phi,method) %>%
                  do({ function(y){
                                  tmp <-  lapply(y, function(x){
                                                     param.family(data.table::fread(x),variable.family)
                                                 })
                                  tmp <- lapply(tmp,function(z){

                                                    z[ceiling(seq(from = 1, to =nrow(z), length.out = length.out )),]
                                                 }
                                  )
                                  as.data.frame(t(colMeans(as.matrix(do.call(rbind,tmp)))))

                                }
                               } (.$directory)
                  )

    }

thin.and.summ(chains,"sigma|mu",500)

#we compute the variance explained
variance.explained <-  function(mcmc.chains, length.out=1001){
         mcmc.chains %>%
             group_by(type,phi,method) %>%
                  do({ function(y){
                                  tmp <-  lapply(y, function(x){
                                                     param.family(data.table::fread(x),"sigma")
                                                 })
                                  tmp <- lapply(tmp,function(z){

                                                    z[ceiling(seq(from = 1, to =nrow(z), length.out = length.out )),]
                                                 }
                                  )
                                  tmp <- lapply(tmp, function(w) {
                                                 w[,'sigmaG[1]']/(rowSums(w))
                                  })
                                  summary.tmp <-  do.call(rbind,tmp)
                                 data.frame(Min  = min(summary.tmp), "0.25.q" = quantile(summary.tmp,0.25), Median = median(summary.tmp), Mean = mean(summary.tmp), "0.75.q" =quantile(summary.tmp,0.75), Max = max(summary.tmp),row.names = NULL)
                                }
                       } (.$directory)
                  )

    }

variance.explained(chains,500)

#we compute the number of markers in the model

markers.in.model <- function(mcmc.chains,length.out){
     mcmc.chains %>%
             group_by(type,phi,method) %>%
                  do({ function(y){
                                  tmp <-  lapply(y, function(x){
                                                     param.family(data.table::fread(x),"comp")
                                                 })
                                   tmp <- lapply(tmp,function(z){

                                                    z[ceiling(seq(from = 1, to =nrow(z), length.out = length.out )),]
                                                 })
                                  tmp <- lapply(tmp,function(z){
                                                     ifelse(z > 0, 1, 0)
                                                 }
                                  )
                                  #data.frame(mean.m.i.m =mean(rowSums(do.call(rbind,tmp))))

                                  summary.tmp <- rowSums(do.call(rbind,tmp))
                                  data.frame(Min  = min(summary.tmp), "0.25.q" = quantile(summary.tmp,0.25), Median = median(summary.tmp), Mean = mean(summary.tmp), "0.75.q" =quantile(summary.tmp,0.75), Max = max(summary.tmp),row.names = NULL)
                                }
                       } (.$directory)
                  )



}
markers.in.model(chains, 500)

pip <- function(mcmc.chains, length.out){
     mcmc.chains %>%
             group_by(type,phi,method) %>%
                  do({ function(y){
                                  tmp <-  lapply(y, function(x){
                                                     param.family(data.table::fread(x),"comp")
                                                 })
                                   tmp <- lapply(tmp,function(z){

                                                    z[ceiling(seq(from = 1, to =nrow(z), length.out = length.out )),]
                                                 })
                                  tmp <- lapply(tmp,function(z){
                                                     ifelse(z > 0, 1, 0)
                                                 }
                                  )
                                  #data.frame(mean.m.i.m =mean(rowSums(do.call(rbind,tmp))))
                                  tmp <- do.call(rbind,tmp)
                                  as.data.frame(t(colSums(tmp)/nrow(tmp)))

                                }
                       } (.$directory)
                  )



}
pip(chains,500)

true.values <-  function() {
   linear.models <-  list.files(data.path,pattern="_linear",full.names =T)
   model.list <- data.frame(linear.models= sapply(linear.models,FUN =basename)) %>% separate("linear.models",sep= "_", into=c("matrix","linear","model"))
   model.list$linear <- NULL
   model.list$model <- NULL
   model.list$directory <- rownames(model.list)
   rownames(model.list) <- NULL

   model.types <- as_tibble(model.list) %>%
           separate( "matrix", into = c("matrix", "phi"), sep = "phi") %>%
           separate( "matrix", into = c("type","variables"), sep ="1k") %>%
           mutate_at(vars(directory),as.character)
   as_tibble(cbind( model.types ,do.call(rbind,lapply(linear.models,function(x){
                                                load(x)
                                                beta.colnames <- sprintf("b[%s]",1:length(b))
                                                tmp.b <- t(b)
                                                colnames(tmp.b) <- beta.colnames
                                                cbind (data.frame(var.eps = var.e, var.g = var.g, var.y =var.y, VE = var.g/var.y), as.data.frame(tmp.b))
                                               }
                                             )
                                )
                )
  )
}

#table with the posterior VE

VE_summary <-  function(mcmc.chains){
     true.values() %>% select(type,phi,VE) %>% inner_join(variance.explained(mcmc.chains,500),by=c("type","phi"))
}
VE_summary(chains)


#now comes the tricky part of putting either the inferred values and the original values in the same scale
true.betas <- true.values() %>% select(type, phi, contains('b['))  %>% gather("coefficient","value", -type, -phi) %>%
  mutate(coefficient = str_extract(coefficient,"[0-9]{1,5}"))



chains.betas <- thin.and.summ(chains,'beta',500)  %>% gather("coefficient","estimate", -type, -phi, -method ) %>%
    mutate(coefficient = str_extract(coefficient, "[0-9]{1,5}"))

#here we compute the RMSE for the betas
true.betas %>% inner_join(chains.betas) %>% group_by(type,phi,method ) %>%  summarise(RMSE = sqrt(mean((value-estimate)^2 )))


#here we do the TP, FP , TN, FN


pip.summ <- function(true.b, chains.b, threshold){

      true.b %>% inner_join(chains.b) %>% group_by(type, phi, method) %>% summarise(TP = sum(ifelse(pip >= threshold & value ==1,1,0 )),
                                                                                                      FP= sum(ifelse(pip >= threshold & value ==0,1,0)),
                                                                                                      TN = sum(ifelse(pip <= threshold & value ==0,1,0)),
                                                                                                      FN = sum(ifelse(pip <= threshold & value ==1,1,0))) %>% mutate(thres = threshold)
}
auc <- function(x,y){
       from <- min(x, na.rm=TRUE)
       to <- max(x, na.rm=TRUE)
       values <- approx(x, y, xout = sort(unique(c(from, to, x[x > from & x < to]))))
       res <- 0.5 * sum(diff(values$x) * (values$y[-1] + values$y[-length(values$y)]))
       res
}
PRcurve<- function(true.b, chains.b){
    thresholds <- seq(0,1,by = 0.01)
    pr <- bind_rows(lapply(thresholds, function(x)pip.summ(true.b, chains.b, x)))
    pr <- pr %>% mutate(precision = TP/(TP+FP), recall = TP/(TP+FN)) %>% unite( "experiment",c(type,phi)) %>% na.omit
    p <- ggplot()+ geom_line(data = pr,aes(x = recall, y= precision, color = experiment, linetype = method )) + theme_bw()+ scale_colour_grey(start = 0, end = .8)
    list(plot = p, table = pr, auc = auc(pr$recall, pr$precision))
}

get.PRcurve <-  function(true.b, chains.pip) {
    true.betas.bool <- true.b %>% mutate(value = ifelse(value !=0 , 1,0))
    chains.betas.bool <- chains.pip %>% gather("coefficient","pip",-type, -phi, -method) %>% mutate(coefficient = str_extract(coefficient,"[0-9]{1,5}"))
    PRcurve(true.betas.bool,chains.betas.bool)
}

tmp <- get.PRcurve(true.betas, pip(chains,500))


#here we are going to plot the  Gelman Rhat for the betas

tmp <- chains.gelman("beta") #note, very expensive operation

gather.tmp<- grep("beta",names(tmp),value=T)


tmp2 <- tmp %>% gather(coefficient,Rhat, gather.tmp) %>% unite("experiment",c(type,phi,method))

#tmp2 <- tmp2 %>% filter(experiment ==  "RAR_0.5_async")

p <-  ggplot() +  geom_histogram(tmp2,mapping =aes(x=Rhat))  +facet_grid(rows=vars(experiment))+ theme_bw()
ggsave("testplot2.pdf",p)

system("/bin/bash extract_running_times.sh")
sync.times <- data.table::fread("sync_times.txt",sep=":")
async.times <- data.table::fread("async_times.txt",sep=":")

colnames(sync.times) <- c("names","hours","minutes","seconds")
colnames(async.times) <- c("names","hours","minutes","seconds")

sync.times <- sync.times %>% transmute(exec_time = 360*hours + 60*minutes + seconds) %>% mutate(method = "sync")
async.times <- async.times %>% transmute(exec_time = 360*hours + 60*minutes + seconds) %>% mutate(method = "async")

p <- ggplot()  + geom_histogram(bind_rows(sync.times,async.times), mapping = aes(x = exec_time),fill = "white" , colour ="black") + theme_bw() + scale_color_grey() + facet_grid(method ~ .)
ggsave("testplot3.pdf",p)


### here we work with the short chains to plot the convergence time
short.files <- list.files(data.path,pattern = "short",full.names = T)
chains.short <- data.frame(directory = short.files, matrices = matrices) %>%
           separate( "matrices", into = c("matrix","method","chain"), sep="_") %>%
           separate( "matrix", into = c("matrix", "phi"), sep = "phi") %>%
           separate( "matrix", into = c("type","variables"), sep ="1k") %>%
           mutate_at(vars(directory),as.character)


variance.explained.trace <-  function(length.out=1001){
         chains.short %>%
             group_by(type,phi,method) %>%
                  do({ function(y){
                                  tmp <-  lapply(y, function(x){
                                                     param.family(data.table::fread(x),"sigma")
                                                 })
                                   #no need to thin for short chains
                                   #tmp <- lapply(tmp,function(z){

                                   #                 z[ceiling(seq(from = 1, to =nrow(z), length.out = length.out )),]
                                    #             }
                                  #)
                                  tmp <- lapply(tmp, function(w) {
                                                 w[,'sigmaG[1]']/(rowSums(w))
                                  })
                                  tmp <- do.call(rbind,tmp)
                                  n.chains <-  nrow(tmp)
                                  colnames(tmp) <- 1:ncol(tmp)
                                  as_tibble(tmp) %>% mutate(chains= 1:n.chains) %>% gather("iteration","VE", -chains)

                                }
                       } (.$directory)
                  )
    }

VE.trace <- variance.explained.trace()
VE.trace <- VE.trace %>% unite(type,c(type,phi)) %>% type_convert(cols(iteration ="i"))
 p <- VE.trace %>% filter(type == "RAR_0.8") %>% ggplot() + geom_path(aes(x=iteration,y=VE,group =chains)) + facet_grid(rows = vars(chains),cols =vars(method)) + theme_bw()
ggsave("testplot4.pdf",p)

trace.async <- VE.trace %>% filter(type == "RAR_0.8", method == "async") %>% group_by(iteration)%>% summarise(mean_VE = mean(VE)) %>% mutate(method = "async")

trace.sync <-  VE.trace %>% filter(type == "RAR_0.8", method == "sync") %>% group_by(iteration)%>% summarise(mean_VE = mean(VE)) %>% mutate(method = "sync")

p <- bind_rows(trace.async, trace.sync) %>% ggplot() + geom_path(aes(x=iteration,y=mean_VE)) + facet_grid(rows =vars(method)) + theme_bw()
ggsave("testplot5.pdf",p)

#plot times for short runs

sync.times <- data.table::fread("sync_times_short.txt",sep=":")
async.times <- data.table::fread("async_times_short.txt",sep=":")

colnames(sync.times) <- c("names","hours","minutes","seconds")
colnames(async.times) <- c("names","hours","minutes","seconds")

sync.times <- sync.times %>% transmute(exec_time = 360*hours + 60*minutes + seconds) %>% mutate(method = "sync")
async.times <- async.times %>% transmute(exec_time = 360*hours + 60*minutes + seconds) %>% mutate(method = "async")

p <- ggplot()  + geom_histogram(bind_rows(sync.times,async.times), mapping = aes(x = exec_time),fill = "white" , colour ="black") + theme_bw() + scale_color_grey() + facet_grid(method ~ .)
ggsave("testplot6.pdf",p)


###############################################################################
##
##               SPEED FILES CODE
##############################################################################


speed.files.path <- "/scratch/temporary/dtrejoba/AsynchronousGibbs"
speed.files <- list.files(speed.files.path,pattern ="speed",full.names =T)

figure.path <- "./PenReg_DiscPrior_Resid_Async/figures/"

get.times  <- function(x){
   sync <- system(paste( c("grep typeppbayes -A1022", x),collapse =" "),intern =T)
   fils  <- sync
   async <- system(paste( c("grep typeasyncppbayes -A1022", x),collapse =" "),intern =T)
   sync <- grep( "Start reading |iteration", sync,value=T)
   async <-  grep( "Start reading |iteration", async,value=T)

   sync <- grep("iteration",sync,value=T)
   async <- grep("iteration",async,value=T)
   sync <- lapply(sync,function(x){strsplit(x,split= ":")[[1]][2] })
   async <- lapply(async,function(x){strsplit(x,split= ":")[[1]][2]})
   sync <- lapply(sync,function(x){
       k <- str_extract_all(x, "\\([^()]+\\)")[[1]] # Get the parenthesis and what is inside
       k <- substring(k, 2, nchar(k)-2)
   })
   async <- lapply(async,function(x){
       k <- str_extract_all(x, "\\([^()]+\\)")[[1]] # Get the parenthesis and what is inside
       k <- substring(k, 2, nchar(k)-2)
   })
   covar <- unlist(lapply(grep("Start reading",fils,value=T),function(x) {rep(strsplit(basename(x),split='X')[[1]][1],1000)}))
   pheno <- unlist(lapply(grep("mcmc-samples",fils, value=T),function(x) {rep(str_extract(x , "_[1-3]s_" )  ,1000) }))
   if(!is.null(covar))
       tibble(iter.times_sync = as.numeric(unlist(sync)),
              iter.times_async = as.numeric(unlist(async)),
              experiment= basename(x) ,
              covariance = covar,
              phenotype = pheno)
   else
      NA
}

execution.times <- lapply(speed.files,get.times)
execution.times[[length(execution.times)]]<- NULL
all.times <- do.call(rbind,execution.times)

all.times <- all.times %>% mutate(covariance = str_sub(covariance,1,4))
all.times <- all.times %>% mutate(covariance = ifelse(covariance == "RAR1", "RAR",covariance))
all.times <- all.times %>% mutate(covariance = ifelse(covariance == "RAR2", "RAR",covariance))
all.times <- all.times %>% separate(experiment,into = c("size", "parallelism"), sep = "speed" )

all.times <- all.times %>% na.omit()


p <- ggplot() + geom_histogram(all.times,mapping=aes(x=iter.times_sync, fill = experiment, alpha =0.2 ))
q <-  ggplot() + geom_histogram(all.times,mapping=aes(x=iter.times_async, fill = experiment, alpha =0.2 ))
library(gridExtra)

ggsave("test_times_sync.pdf",p)
ggsave("test_times_async.pdf",q)

r <- ggplot() + geom_histogram(all.times, mapping =aes(x=iter.times_sync/iter.times_async,fill = experiment, alpha =0.2))
ggsave("test_times_speedup.pdf",r)

#this plot is cool but it may be better to plot mean iteration times.
s <- all.times %>% ggplot(aes(x=iter.times_sync / iter.times_async, fill = parallelism)) + geom_histogram(aes(y=0.5*..density..),binwidth=0.5) + facet_grid(cols = vars(size), rows=vars(covariance))+geom_vline(xintercept = 1)
ggsave("test_times_speedup_par.pdf",s)



 mean.times <- all.times %>% group_by(size, parallelism,covariance, phenotype) %>% summarise(mean_sync= mean(iter.times_sync), mean_async = mean(iter.times_async))

speed.ups <- mean.times %>% mutate(speed_up = mean_sync / mean_async)

t.p <-  speed.ups %>% ggplot(aes(x= covariance, y = speed_up ,colour = parallelism)) + geom_point() + facet_grid( rows = vars(size) )
ggsave("test_mean_speedup_par.pdf",t.p)

single.threaded <- all.times %>% group_by(size,covariance,phenotype) %>%  summarise(single.thread = mean(iter.times_sync[parallelism =="_1_1_1_1"]))

single.speedup <- mean.times %>%
    inner_join(single.threaded) %>%
    mutate(sync = ( single.thread /mean_sync) ,async = (single.thread/mean_async )) %>%
    gather(type, speed.up, c(sync,async))


#this is the good plot !
t.q <- single.speedup %>%
    ggplot(aes(x = covariance, y = speed.up, colour =parallelism, shape =type)) +
    geom_point() +
    facet_grid(rows =vars(size),cols = vars(phenotype)) +
    scale_shape_manual(c(3,4))

ggsave(paste0(figure.path,"single_speedup_par.pdf"),t.q)
