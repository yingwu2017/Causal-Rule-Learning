library(grf)
library(dplyr)
library(pre)
library(glmnet)
library(ggplot2)

path <- 'simulation/study two/'
data.path <- paste0(path,'data/')
source('functions/rule.process.R', encoding = 'UTF-8')
source('functions/CRL_functions.R', encoding = 'UTF-8')


# define settings
Ns <- c(2500,5000,10000) 
eff.sizes <- seq(5,30,5)
num.trs <- c(10,30,50)
num.crs <- c(100,300,500) # number of candidate rules

# result container
d.rule.codes.ls <- list()
for (num.cr in num.crs) {
  for (N in Ns) {
    for (eff.size in eff.sizes) {
      for (num.tr in num.trs) {
        
        file <- paste0(N,'obs;',eff.size,'tau;',num.cr,'num.candi;',num.tr,'true_grp')
        # load data
        datafile <- paste0(N,'obs;',eff.size,'tau;',num.cr,'num.candi;',num.tr,'true_grp','.RData')
        paste0('This is for data: ',datafile) %>% print()
        load(paste0(data.path,datafile))
        
        # define propensity scores
        ps <- rep(0.5,nrow(data))
        
        rule.bi <- data[1:num.cr]
        
        # ---- train D-learning model ----
        y.new <- 2*data$outcome*data$treat
        X <- as.matrix(rule.bi) 
        weights <- 1/ps
        # tune parameters
        cv.rule.model <- cv.glmnet(x = X,
                                   y = y.new,
                                   weights = weights,
                                   type.measure = 'mse',
                                   nfolds = 10,
                                   trace.it = 1,
                                   seed= 2025)
        
        
        if(is.null(cv.rule.model)) next()
        d.rule.model <- glmnet(x = X,
                               y = y.new,
                               family = 'gaussian',
                               weights = weights, 
                               alpha = 1, 
                               lambda = cv.rule.model$lambda.min,
                               standardize = F,
                               intercept = F,
                               trace.it = T)
        d.rule.codes.ls[[file]] <- colnames(X)[which(d.rule.model$beta!=0)] %>% gsub(pattern = 'rule',replacement = '') %>% as.integer()
      }
    }
  }
} 



#------ visualization of results ------
#------ construct data frame to plot: Number of non-zero rules selected by D-learning.----
df <- data.frame()
for (num.cr in num.crs){
  for (N in Ns){
    for (eff.size in eff.sizes){
      index <- paste0(N,'obs;',eff.size,'tau;',num.cr,'num.candi;',num.trs,'true_grp')
      values <- lapply(index, FUN = function(ind){
        length(d.rule.codes.ls[[ind]])
      }) %>% unlist()
      temp <- data.frame(M = rep(num.cr,3)
                         ,k = rep(eff.size,3)
                         ,N = rep(N,3)
                         ,num.true.grp = num.trs
                         ,num.nzero.rules = values
      )
      df = rbind(df,temp)
      
    }
  }
}

df$num.true.grp <- as.factor(df$num.true.grp)
df$N <- as.factor(df$N)
df$M <- as.factor(df$M)
df$k <- as.factor(df$k)

color.values <- c("10" = 'gray', '30' = '#FC7C91','50' = '#218db1')
ggplot(df, aes(x = N, y = num.nzero.rules, fill = num.true.grp)) +
  geom_col(width = .5, 
           position = position_dodge(width = .5)) +  
  geom_text(aes(label = num.nzero.rules)
            ,position = position_dodge(width = .5)
            ,vjust = 2
            ,size = 2                              
            ,color = "black") +  
  scale_fill_manual(values = color.values)+
  facet_grid(k~M,labeller = label_both)+
  labs(x =  "N", y = 'Number of Non-zero Rules selected by D-learning')+  
  theme_bw()+  
  theme(legend.position = "none")

#------ construct data frame to plot: Number of True Rules Identified by D-learning.----
df2 <- data.frame()
for (num.tr in num.trs){
  true.inds <- if (num.tr == 10) {
    1:10
  } else if (num.tr == 30) {
    1:30
  } else if (num.tr == 50) {
    1:50
  }
  for (num.cr in num.crs){
    for (N in Ns){
      for (eff.size in eff.sizes){
        index <- paste0(N,'obs;',eff.size,'tau;',num.cr,'num.candi;',num.tr,'true_grp')
        rule.inds <- lapply(index, FUN = function(ind){
          d.rule.codes.ls[[ind]]
        }) %>% unlist()
        value <- length(intersect(rule.inds,true.inds))
        with(df, which(M == num.cr & k == eff.size & num.true.grp == num.tr & N == N))
        
        temp <- data.frame(M = num.cr
                           ,k = eff.size
                           ,N = N
                           ,num.true.grp = num.tr
                           ,num.true.found = value
        )
        df2 = rbind(df2,temp)
      }
    }
  }
}


df2$num.true.grp <- as.factor(df2$num.true.grp)
df2$N <- as.factor(df2$N)
df2$M <- as.factor(df2$M)
df2$k <- as.factor(df2$k)

color.values <- c("10" = 'gray', '30' = '#FC7C91','50' = '#6AA59D')
ggplot(df2, aes(x = N, y = num.true.found, fill = num.true.grp)) +
  geom_col(width = .5, 
           position = position_dodge(width = .5)) +  
  geom_text(aes(label = num.true.found)
            ,position = position_dodge(width = .5)
            ,vjust = 2
            ,size = 2                              
            ,color = "black") +  
  scale_fill_manual(values = color.values)+
  facet_grid(k~M,labeller = label_both)+
  labs(x =  "N", y = 'Number of True Rules Identified by D-learning')+  
  theme_bw()+  
  theme(legend.position = "none")


