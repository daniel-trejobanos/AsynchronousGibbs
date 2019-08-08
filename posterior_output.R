## TODO : Add the level of parallelism as a grouping variable and add this facet to the plots
##' Posterior analysis code
##' First we select the data path, lets start with a single data set size
data.path <- "/scratch/temporary/dtrejoba/AsynchronousGibbs/data/2k5k"

##' We read the short runs only
files.path <- list.files(data.path, pattern ="s_.*_short", full.names = T)
matrices <-   list.files(data.path, pattern = "s_.*_short")

require(coda)
require(tidyverse)
##' We still havent had the time to properly implement simulation based callibration
#source('./sbc.R')

##' We save ta data frame with the file paths and descriptors of the experiment
chains <- data.frame(directory = files.path, matrices = matrices) %>%
    separate( "matrices", into = c("matrix","sparsity","method","chain","chain-type"), sep="_") %>%
    separate( "matrix", into = c("matrix", "phi"), sep = "phi") %>%
    separate( "matrix", into = c("type","variables"), sep ="1k") %>%
    mutate_at(vars(directory),as.character) %>%
    as_tibble()

##' function to extract a family of parameters from a mcmc chains
param.family <- function( chain.matrix, family){
#  if(is.character(family))
   as.mcmc(select(chain.matrix,matches(family)))
 # else
  #  stop("error, family has to be character")
}

##' we create a path to save the figures or tempoary files from the posterior summary pipeline
posterior.path <- paste(data.path,"/post",sep ="")
dir.create(posterior.path)

##' Function to estimate the gelman R for each of the parameters in the parameter family
chains.gelman <- function(mcmc.chains,variable.family, multivariate = F){
                  mcmc.chains %>%
                  group_by(type,phi,method) %>%
                  do({ function(y){
                                  tmp <- coda::gelman.diag( mcmc.list(lapply(y, function(x){
                                                     param.family(data.table::fread(x),variable.family)
                                                 })),
                                                autoburnin =T,
                                                multivariate = multivariate)
                                   tmp.df <- as.data.frame(tmp$psrf)
                                   tmp.df$var <- rownames(tmp$psrf)
                                   tmp.df$"Upper C.I." <- NULL
                                   tmp.df %>% spread(var, "Point est.")

                                }
                               } (.$directory)
                  )



                 }


##' we compute the effective sample size for each
#tmp<-chains.gelman(chains,"sigma|mu")

##' function to compute the effective sample size for the chains in the mcmcm.chains tibble for all the parameters of the given family
chains.effective <-  function(mcmc.chains,variable.family, burnin){
    tmp <- lapply(mcmc.chains$directory, FUN = function(x) { coda::effectiveSize(param.family(data.table::fread(x),
                                                                                              variable.family)[burnin:3999, ])})
                  as_tibble(cbind(chains,do.call(rbind,tmp)))
                 }

##' function to compute the average effective sample size over chains in the same experiment
avg.effective <- function(mcmc.chains, variable.family,burnin){
                  chains.effective(mcmc.chains, variable.family,burnin) %>%
                     group_by(type,phi,method,sparsity) %>%
                      summarise_if(is.numeric , c(mean,sd))
    }

#tmp <-  avg.effective(chains, "sigma|mu",2000)

##' Here we want to read the chains and thin them and compute the parameters summaries
thin.and.summ <- function(mcmc.chains,variable.family,length.out=1001,burnin){
         mcmc.chains %>%
             group_by(type,phi,method,sparsity) %>%
                  do({ function(y){
                                  tmp <-  lapply(y, function(x){
                                                     param.family(data.table::fread(x),variable.family)[burnin:3999,]
                                                 })
                                  tmp <- lapply(tmp,function(z){

                                                    z[ceiling(seq(from = 1, to =nrow(z), length.out = length.out )),]
                                                 }
                                  )
                                  as_tibble(t(colMeans(as.matrix(do.call(rbind,tmp)))))

                                }
                               } (.$directory)
                  )

    }

#thin.and.summ(chains,"sigma|mu",500,2000)

##' we compute the variance explained and its statistics for all the chains in the data frame
variance.explained <-  function(mcmc.chains, length.out=1001,burnin){
         mcmc.chains %>%
             group_by(type,phi,method,sparsity) %>%
                  do({ function(y){
                                  tmp <-  lapply(y, function(x){
                                                     param.family(data.table::fread(x),"sigma")[burnin:3999,]
                                                 })
                                  Tmp <- lapply(tmp,function(z){

                                                    z[ceiling(seq(from = 1, to =nrow(z), length.out = length.out )),]
                                                 }
                                  )
                                  tmp <- lapply(tmp, function(w) {
                                                 w[,'sigmaG[1]']/(rowSums(w))
                                  })
                                  summary.tmp <-  do.call(rbind,tmp)
                                 tibble(Min  = min(summary.tmp), "0.25.q" = quantile(summary.tmp,0.25), Median = median(summary.tmp), Mean = mean(summary.tmp), "0.75.q" =quantile(summary.tmp,0.75), Max = max(summary.tmp))
                                }
                       } (.$directory)
                  )

}

#variance.explained(chains,500,2000)

##' we compute the number of markers in the model
markers.in.model <- function(mcmc.chains,length.out,burnin){
     mcmc.chains %>%
             group_by(type,phi,method,sparsity) %>%
                  do({ function(y){
                                  tmp <-  lapply(y, function(x){
                                                     param.family(data.table::fread(x),"comp")[burnin:3999,]
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
                                  tibble(Min  = min(summary.tmp), "0.25.q" = quantile(summary.tmp,0.25), Median = median(summary.tmp), Mean = mean(summary.tmp), "0.75.q" =quantile(summary.tmp,0.75), Max = max(summary.tmp))
                                }
                       } (.$directory)
                  )



}
#markers.in.model(chains, 500,2000)

##' this functions yields a tibble with markers as columns and the PIP as values
pip <- function(mcmc.chains, length.out,burnin){
     mcmc.chains %>%
             group_by(type,phi,method,sparsity) %>%
                  do({ function(y){
                                  tmp <-  lapply(y, function(x){
                                                     param.family(data.table::fread(x),"comp")[burnin:3999,]
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
                                  as_tibble(data.frame(t(colSums(tmp)/nrow(tmp))))

                                }
                       } (.$directory)
                  )



}
#pip(chains,500,2000)

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
                                                beta.colnames <- sprintf("b[%s]",1:ncol(b))
                                                tmp.b <- b
                                                colnames(tmp.b) <- beta.colnames
                                                as_tibble(cbind (data.frame(sparsity =c("1s","2s","3s") ,var.eps = var.e, var.g = var.g, var.y =var.y, VE = var.g/var.y), as.data.frame(tmp.b)))
                                               }
                                             )
                                )
                )
  )
}

##' table with the ve and the summary statistics of the thinned and post burn in samples
VE_summary <-  function(mcmc.chains){
     true.values() %>% select(type,phi,sparsity,VE) %>% inner_join(variance.explained(mcmc.chains,500,2000),By=c("type","phi","sparsity"))
}
#VE_summary(chains)


##' this dplyr query gets the true betas in tall format
true.betas <- true.values() %>% select(type, phi, sparsity, var.y ,contains('b['))  %>% gather("coefficient","value", -type, -phi , -sparsity, -var.y) %>%
  mutate(coefficient = str_extract(coefficient,"[0-9]{1,5}"))



chains.betas <- thin.and.summ(chains,'beta',500,2000)  %>% gather("coefficient","estimate", -type, -phi, -method, -sparsity ) %>%
    mutate(coefficient = str_extract(coefficient, "[0-9]{1,5}"))

#' here we compute the RMSE for the betas
true.betas %>% inner_join(chains.betas) %>%
    group_by(type,phi,method,sparsity) %>%
    summarise(RMSE = sqrt(mean((value-var.y*estimate)^2 ))) %>%


#here we do the TP, FP , TN, FN


pip.summ <- function(true.b, chains.b, threshold){

    true.b %>% inner_join(chains.b) %>%
        group_by(type, phi, method,sparsity) %>%
        summarise(TP = sum(ifelse(pip >= threshold & value ==1,1,0 )),
                  FP= sum(ifelse(pip >= threshold & value ==0,1,0)),
                  TN = sum(ifelse(pip <= threshold & value ==0,1,0)),
                  FN = sum(ifelse(pip <= threshold & value ==1,1,0))) %>%
          mutate(thres = threshold)
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
    pr <- pr %>% mutate(precision = TP/(TP+FP), recall = TP/(TP+FN), no_skill = ifelse(sparsity =="1s",0.1,ifelse(sparsity=="2s",0.5,0.9)) ) %>%
        unite( "experiment",c(type,phi)) %>% na.omit()
    p <- pr %>%
        ggplot()+
        geom_line(aes(x = recall, y= precision, color = experiment, linetype = method )) +
        theme_bw()+ scale_colour_grey() +
        facet_grid(rows = vars(sparsity) )+
        geom_hline(aes(yintercept = no_skill))


    list(plot = p, table = pr, auc = auc(pr$recall, pr$precision))
}

get.PRcurve <-  function(true.b, chains.pip) {
    true.betas.bool <- true.b %>% mutate(value = ifelse(value !=0 , 1,0))
    chains.betas.bool <- chains.pip %>% gather("coefficient","pip",-type, -phi, -method, -sparsity) %>% mutate(coefficient = str_extract(coefficient,"[0-9]{1,5}"))
    PRcurve(true.betas.bool,chains.betas.bool)
}

##' we plot the PR curves faceted by the sparsity level
tmp <- get.PRcurve(true.betas, pip(chains,500,2000))


#here we are going to plot the  Gelman Rhat for the betas

tmp <- chains.gelman("beta") #note, very expensive operation

gather.tmp<- grep("beta",names(tmp),value=T)


tmp2 <- tmp %>% gather(coefficient,Rhat, gather.tmp) %>% unite("experiment",c(type,phi,method))

#tmp2 <- tmp2 %>% filter(experiment ==  "RAR_0.5_async")

p <-  ggplot() +  geom_histogram(tmp2,mapping =aes(x=Rhat))  +facet_grid(rows=vars(experiment))+ theme_bw()
ggsave("testplot2.pdf",p)




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



mean.times <- all.times %>% group_by(size, parallelism,covariance, phenotype) %>%
    summarise(mean_sync= mean(iter.times_sync), mean_async = mean(iter.times_async))

speed.ups <- mean.times %>%
    mutate(speed_up = mean_sync / mean_async)

single.threaded <- all.times %>%
    group_by(size,covariance,phenotype) %>%
    summarise(single.thread = mean(iter.times_sync[parallelism =="_1_1_1_1"]))

single.speedup <- mean.times %>%
    inner_join(single.threaded) %>%
    mutate(sync = ( single.thread /mean_sync) ,async = (single.thread/mean_async )) %>%
    gather(type, speed.up, c(sync,async)) %>%
    ungroup() %>%
    mutate( parallelism = fct_recode(parallelism,
                                 "1" = "_1_1_1_1",
                                 "2" = "_2_2_2_2",
                                 "4" = "_4_4_4_4",
                                 "8" = "_8_8_8_8",
                                 "16" = "_16_16_16_16",
                                 "32" = "_32_32_32_32")) %>%
    mutate(parallelism = as.numeric(as.character(parallelism)))%>%
    mutate(size = fct_recode(size,
                             "1"="1k1k_",
                             "2"="1k5k_",
                             "3"="2k5k_",
                             "4"="10k50k_")) %>%
    mutate_at(vars(size), funs(factor(., levels=unique(.))))

single.speedup$size <- factor(single.speedup$size, levels = c("1","2","3","4"))


#this is the good plot !

s.labels <- c("10% non-zero", "50% non-zero", "90% non-zero")
names(s.labels) <- c("_1s_","_2s_","_3s_")

size.labels <- c("1e3x1e3","1e3x5e3","2e3x5e3","1e4x5e4")
names(size.labels) <- c("1","2","3","4")


t.q  <-single.speedup %>%
    ggplot(aes(x = parallelism, y = speed.up, colour =covariance, shape =type)) +
      geom_hline(yintercept = 1) +
     scale_shape_manual(values= c(3,4)) +
    geom_line(aes(linetype = type)) +
    scale_color_grey() +
    facet_grid(rows =vars(size),cols = vars(phenotype), labeller= labeller(phenotype=s.labels, size = size.labels))+
    theme_bw() +
     theme(panel.border = element_blank(), panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    xlab("Threads")+
    ylab(expression("Speedup(",frac("m.i.t. parallel"," m.i.t. single thread"),")")))







ggsave(paste0(figure.path,"single_speedup_par.pdf"),t.q)


