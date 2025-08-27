
#-------load required packages---
library(dplyr)
#-------define data settings-------

# define the varying sample sizes
Ns <- c(2500,5000,10000)

# define the varying number of true subgroups
num.true.grps <- c(1,3,5)

# define the varying number of effect size quantity k
eff.sizes <- seq(5,30,5)

# define the varying number of covariates
num.covar <- 9  

# number of resampled datasets for each data setting
num.run <- 100 

# generate simulation data
for (N in Ns){
  for (num.true.grp in num.true.grps) {
    for (eff.size in eff.sizes){
      #-------data generation for a given setting-------
      X <- data.frame(V1 = rbinom(n = N, size = 1, prob = 0.5),
                      V2 = rbinom(n = N, size = 1, prob = 0.5))
      if(num.true.grp == 5){
        X[,3] <- runif(n = N, min = 0, max = 10)
      } else {
        X[,3] <- rbinom(n = N, size = 1, prob = 0.5)
      }
      for(i in 4:num.covar){ 
        X[,i] <- rnorm(n = N, mean = 0, sd = 2)
      }
      colnames(X) <- paste0('x',1:num.covar) 
      
      # treatment assignment
      treat <- rbinom(n = N, size = 1, prob = 0.5)
      treat[which(treat==0)] <- -1
      X$treat <- treat
      
      
      attach(X)
      # define main effects:
      main.effect <- 1+2*x1**2+3*x1*x2+0.4*x4*x7+x5-x6+0.5*x3*x8
      
      # define true subgroups and corresponding treatment effects
      if(T){
        if(num.true.grp == 1){
          true.cate <- rep(0, N)
          true.cate[which(x1==0 & x2==0)] <- -eff.size
          true.cate[which(x1==1 & x2==1)] <- eff.size
          summary(true.cate)
        }
        if(num.true.grp == 3){
          true.cate <- rep(0, N)
          true.cate[which(x1==1 & x2==0 & x3==0)] <- eff.size
          true.cate[which(x1==1 & x2==1 & x3==0)] <- 2*eff.size
          true.cate[which(x1==1 & x2==1 & x3==1)] <- 3*eff.size
          true.cate[which(x1==0 & x2==0 & x2==0)] <- -eff.size
          summary(true.cate)
        }
        if(num.true.grp == 5){
          true.cate <- rep(0, N)
          true.cate[which(x1==1 & x2==1 & x3>5)] <- 3*eff.size
          true.cate[which(x1==0 & x2==1 & x3>5)] <- 2*eff.size
          true.cate[which(x1==1 & x2==0 & x3>5)] <- 2*eff.size
          true.cate[which(x1==1 & x2==1 & x3<=5)] <- 2*eff.size
          true.cate[which(x1==0 & x2==0 & x3>5)] <- eff.size
          true.cate[which(x1==0 & x2==0 & x3<=5)] <- -eff.size
          summary(true.cate)
        }
        
      }   
      
      detach(X)
      
      # generate potential outcomes
      
      # continuous and positive
      y_0 <- main.effect + rnorm(N,mean = 0,sd = 1)
      y_1 <- y_0+true.cate
      
      # assign observed outcome
      X[X$treat==-1,'outcome'] <-  y_0[X$treat==-1]
      X[X$treat== 1,'outcome'] <-  y_1[X$treat== 1]
      
      data <- X
      
      
      #-------re-sample data and then split into training and test data-------
      data.train <- list()
      data.test <- list()
      
      y_0.tr <- list()
      y_0.te <- list()
      
      y_1.tr <- list()
      y_1.te <- list()
      
      true.cate.tr <- list()
      true.cate.te <- list()
      
      
      
      for(i in 1:num.run){
        set.seed(i) 
        index <- sample(0:1,nrow(data),replace = T, prob = c(0.3,0.7))
        
        tr.samp <- which(index == 1)
        te.samp <- which(index == 0)
        
        data.train[[i]] <- data[tr.samp,]
        data.test[[i]]  <- data[te.samp,]
        
        y_0.tr[[i]] <- y_0[tr.samp]
        y_1.tr[[i]] <- y_1[tr.samp]
        
        y_0.te[[i]] <- y_0[te.samp]
        y_1.te[[i]] <- y_1[te.samp]
        
        true.cate.tr[[i]] <- true.cate[tr.samp]
        true.cate.te[[i]] <- true.cate[te.samp]
      }
      
      #-------save all re-sampled training and test datasets-------- 
      save(data.train, data.test, 
           y_0,y_1,y_0.tr,y_1.tr,y_0.te,y_1.te,
           true.cate,true.cate.tr,true.cate.te,
           file = paste0('simulation/study one/data/',N,'obs;',eff.size,'tau;',num.true.grp,'grps','.RData')
      ) 
    } 
  }
}


