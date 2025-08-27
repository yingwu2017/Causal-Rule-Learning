library(dplyr)
path <- 'simulation/study two/'
data.path <- paste0(path,'data/')

#-------define data settings-------
Ns <- c(2500,5000,10000) 
eff.sizes <- seq(5,30,5)
num.trs <- c(10,30,50) # number of true rules
num.crs <- c(100,300,500) # number of candidate rules

set.seed(2022)
#-------data generation-------

# generate covariates
for (num.cr in num.crs) {
  for (N in Ns) {
    for (eff.size in eff.sizes){
      for (num.tr in num.trs){
        X <- matrix(data = rep(0,N*num.cr),nrow = N) %>% as.data.frame()
        
        probs <- runif(num.cr,min = 0.1,max = 0.9)
        for (i in 1:num.cr){
          X[i] <- rbinom(N, size = 1, prob = probs[i])
        }
        colnames(X) <- paste0('rule',1:num.cr)
        # define initial subgroup cate at random
        subgrps.cate <- rnorm(num.cr, mean = 0, sd = 1) #rep(0,num.cr)
        # define the coefficients for the true cate
        cf <- rep(0, num.cr)
        if(num.tr == 10){
          true.index <- 1:10
          cf[true.index] <- .1
        }
        if(num.tr == 30){
          true.index <- 1:30
          cf[true.index] <- 1/30
        }
        
        if(num.tr == 50){
          true.index <- 1:50
          cf[true.index] <- 1/50
        }
        subgrps.cate[true.index] <- eff.size
        # transform the 0-1 matrix to 0-CATEs
        for (i in 1:num.cr){
          X[,i] <- ifelse(X[,i] == 1, subgrps.cate[i],0)
        }
        
        # compute the true cate
        attach(X)
        true.rules <- paste0('rule',true.index)
        order <- paste(cf,colnames(X), sep = '*', collapse = '+')
        true.cate <- parse(text = order) %>% eval()
        detach(X)
        
        # treatment assignment
        treat <- rbinom(n = N, size = 1, prob = 0.5)
        
        table(treat) 
        treat[which(treat==0)] <- -1
        X$treat <- treat
        
        # generate potential outcome
        y_0 <- rnorm(N, mean = 0, sd = eff.size/10)
        y_1 <- y_0+true.cate
        
        # assign observed outcome
        X[X$treat==-1,'outcome'] <-  y_0[X$treat==-1]
        X[X$treat== 1,'outcome'] <-  y_1[X$treat== 1]
        
        data <- X
        
        data.filename <- paste0(N,'obs;',eff.size,'tau;',num.cr,'num.candi;',num.tr,'true_grp','.RData')
        save(data, 
             y_0,y_1,
             true.cate,
             file = paste0(data.path,
                           data.filename)
        )    
        
      }
      
    }
  }  
}






