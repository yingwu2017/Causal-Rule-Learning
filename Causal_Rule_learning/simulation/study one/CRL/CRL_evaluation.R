
library(dplyr)
path <- 'simulation/study one/'

# you must source this .R file before you run the codes below 
source('functions/rule.process.R', encoding = 'UTF-8')


# create containers to save performance metrics result 
cr.MeanPot <- list()
cr.MeanMse <- list()
cr.MeanOverlap <- list()

cr.num.find <- list()
cr.num.fake <- list()
cr.prop.true <- list()


#---- define the settings to be evaluated -----

Ns <- c(2500,5000,10000) 
true.num.grps <- c(1,3,5)
eff.sizes <- seq(5,30,5)
num.run <- 100



#----the evaluation loop begins---- 
for (N in Ns){
  # N <- 2500
  # eff.size <- 10
  # num.true.grp <- 1
  for (eff.size in eff.sizes) {
    for (num.grp in true.num.grps){
      #-------define the true rules---------
      if(num.grp==1){
        true.rules <- c("x1==1 & x2==1"
                        ,"x1==0 & x2==0"
        )
      }
      if(num.grp==3){
        true.rules <- c( "x1==1 & x2==0 & x3==0"
                         ,"x1==1 & x2==1 & x3==0"
                         ,"x1==1 & x2==1 & x3==1" 
                         ,"x1==0 & x2==0 & x3==0"
        )
      }
      if(num.grp==5){
        true.rules <- c( "x1==1 & x2==1 & x3>5"
                         ,"x1==0 & x2==1 & x3>5"
                         ,"x1==1 & x2==0 & x3>5"
                         ,"x1==1 & x2==1 & x3<=5"
                         ,"x1==0 & x2==0 & x3>5"
                         ,"x1==0 & x2==0 & x3<=5"
        )
      }
      
      # load data and fitted model ---
      datafile <- paste0(N,'obs;',eff.size,'tau;',num.grp,'grps.RData')
      load(paste0(path,'data/',datafile)) 
      
      resultfile <- paste0(N,'obs;',eff.size,'tau;',num.grp,'grps_result.RData')
      load(paste0(path,'CRL/model.result/',resultfile))
      
      
      ##---- performance metrics container ----
      mean.PotentialOutcomes <- c()
      MSE.cate <- c()
      overlap.rate <- c()
      num.fake <- c()
      num.true_grp <- c()
      prop.true <- c()
      
      for (i in 1:num.run){
        #i <- 1
        print(i)
        ##----load model results----
        cate.estimate <- cate.test[[i]]
        predictions  <-  ifelse(cate.estimate>0,1,-1)
        estimate.grp <- which(predictions==1)
        
        test <- data.test[[i]]
        num.sample <-nrow(test)
        true.grp <- which(true.cate.te[[i]] > 0)
        
        ##-- compute the metrics----
        #----- MSE, MPO, Overlap ----
        
        #MSE of cate estimates
        MSE.cate[i] <-  compute.Mse.CateEstimates(true.cate.te[[i]],cate.estimate)
        
        #MPO of recommended treatment
        treatA_outcomes <- y_0.te[[i]][which(predictions==-1)]
        treatB_outcomes <- y_1.te[[i]][which(predictions==1)]
        cur.val <- (sum(treatA_outcomes)+sum(treatB_outcomes))/length(predictions)
        mean.PotentialOutcomes[i] <- round(cur.val,digits = 3) #c(mean.PotentialOutcomes,)
        
        # overlap rate
        overlap.rate[i] <- length(intersect(estimate.grp,true.grp))/length(union(estimate.grp,true.grp))
        
        
        
        
        ##---- reformat the representations of identified rules----
        
        
        Subgrps[[i]]$rule <- lapply(Subgrps[[i]]$rule, 
                                    FUN = rule.thre.adjust,
                                    var.name = 'x3',
                                    var.thre = 5, 
                                    acc.ran = 0.5) %>% unlist()
        
        Subgrps[[i]]$rule <- lapply(Subgrps[[i]]$rule, 
                                    FUN = format_CRL_rule) %>% unlist()
        
        Subgrps[[i]] <- rm.duplicates(Subgrps[[i]])
        
        #count the number of true/fake groups identified by each run
        num.true_grp[i] <- lapply(unique(Subgrps[[i]]$rule), FUN = function(rule){
          ifelse(rule %in% true.rules,1,0)}) %>% unlist() %>% sum()
        num.fake[i] <- nrow(Subgrps[[i]])-num.true_grp[i]
        # proportion of the true rules in the identified rule set
        prop.true[i] <- num.true_grp[i]/nrow(Subgrps[[i]])
      }
      
      
      
      
      
      
      ##---- save evaluation results of 100 runs----
      save.index <- paste0(N,'_',eff.size,'_',num.grp)
      
      cr.MeanPot[[save.index]] <- mean.PotentialOutcomes
      cr.MeanMse[[save.index]] <- MSE.cate
      cr.MeanOverlap[[save.index]] <- overlap.rate
      
      cr.num.find[[save.index]] <- num.true_grp
      cr.num.fake[[save.index]] <- num.fake
      cr.prop.true[[save.index]] <- prop.true
      
    }
  }
}

#----save results----
save(cr.MeanMse,
     cr.MeanPot, 
     cr.MeanOverlap, 
     cr.num.find,
     cr.num.fake,
     cr.prop.true,
     file = paste0(path,'CRL/performance_metrics_CRL.RData')
)



