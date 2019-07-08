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
chains.gelman <- function(variable.family, multivariate = F){
                  chains %>%
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
tmp<-chains.gelman("sigma|mu")
chains.effective <-  function(variable.family){
                  tmp <- lapply(chains$directory, FUN = function(x) { coda::effectiveSize(param.family(data.table::fread(x),variable.family))})  
                  cbind(chains,do.call(rbind,tmp))
                 }

chains.effective("sigma|mu")

avg.effective <- function(variable.family){
                  chains.effective(variable.family) %>%
                     group_by(type,phi,method) %>%
                      summarise_if(is.numeric , c(mean,sd))                      
    }

tmp <-  avg.effective("sigma|mu")

# Here we want to read the chains and thin them and compute the parameters summaries
thin.and.summ <- function(variable.family,length.out=1001){
         chains %>%
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

thin.and.summ("sigma|mu",500)

#we compute the variance explained
variance.explained <-  function(length.out=1001){
         chains %>%
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

variance.explained(500)

#we compute the number of markers in the model

markers.in.model <- function(length.out){
     chains %>%
             group_by(type,phi,method) %>%
                  do({ function(y){
                                  tmp <-  lapply(y, function(x){
                                                     param.family(data.table::fread(x),"comp")
                                                 })
                                   tmp <- lapply(tmp,function(z){
                                      
                                                    z[ceiling(seq(from = 1, to =nrow(z), length.out = length.out )),]
                                                 }
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
markers.in.model(500)

pip <- function(length.out){
     chains %>%
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
pip(500)

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
true.values() %>% select(type,phi,VE) %>% inner_join(variance.explained(500),by=c("type","phi"))

#now comes the tricky part of putting either the inferred values and the original values in the same scale
true.betas <- true.values() %>% select(type, phi, contains('b[')) 
gather.col.true <- grep('b',names(true.betas),value=T)

true.betas <- true.betas %>% gather("coefficient","value", gather.col.true) %>% 
  mutate(coefficient = str_extract(coefficient,"[0-9]{1,5}"))



chains.betas <- thin.and.summ('beta',500)

gather.cols.chains <- grep('beta', names(chains.betas), value =T)

chains.betas <- chains.betas %>% gather("coefficient","estimate", gather.cols.chains) %>% mutate(coefficient = str_extract(coefficient, "[0-9]{1,5}"))

#here we compute the RMSE for the betas
true.betas %>% inner_join(chains.betas) %>% group_by(type,phi,method ) %>%  summarise(RMSE = sqrt(mean((value-estimate)^2 )))


#here we do the TP, FP , TN, FN
true.betas.bool <- true.betas %>% mutate(value = ifelse(value !=0 , 1,0))
chains.betas.bool <- pip(500)
gather.cols.chains.bool <- grep('comp', names(chains.betas.bool), value =T)
chains.betas.bool <- chains.betas.bool %>% gather("coefficient","pip",gather.cols.chains.bool) %>% mutate(coefficient = str_extract(coefficient,"[0-9]{1,5}"))

pip.summ <- function(threshold){
      
      true.betas.bool %>% inner_join(chains.betas.bool) %>% group_by(type, phi, method) %>% summarise(TP = sum(ifelse(pip >= threshold & value ==1,1,0 )),fsaf
                                                                                                      FP= sum(ifelse(pip >= threshold & value ==0,1,0)),
                                                                                                      TN = sum(ifelse(pip <= threshold & value ==0,1,0)),
                                                                                                      FN = sum(ifelse(pip <= threshold & value ==1,1,0))) %>% mutate(thres = threshold)
}

PRcurve <- function(){
    thresholds <- seq(0,1,by = 0.01)
    pr <- bind_rows(lapply(thresholds, function(x)pip.summ(x)))
    pr <- pr %>% mutate(precision = TP/(TP+FP), recall = TP/(TP+FN)) %>% unite( "experiment",c(type,phi))
    p <- ggplot()+ geom_line(data = pr,aes(x = recall, y= precision, color = experiment, linetype = method )) + theme_bw()+ scale_colour_grey(start = 0, end = .8)
    list(plot = p, table = pr)
}
     
#here we are going to plot the  Gelman Rhat for the betas

tmp <- chains.gelman("beta") #note, very expensive operation

gather.tmp<- grep("beta",names(tmp),value=T)

     
tmp2 <- tmp %>% gather(coefficient,Rhat, gather.tmp) %>% unite("experiment",c(type,phi,method)) 

#tmp2 <- tmp2 %>% filter(experiment ==  "RAR_0.5_async")

p <-  ggplot()+  geom_histogram(tmp2,mapping =aes(x=Rhat))  +facet_grid(rows=vars(experiment))+ theme_bw()
ggsave("testplot2.pdf",p)
    
