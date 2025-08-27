library(grf)
library(dplyr)
library(glmnet)

# you must source this .R file before you run CRL models  
source(paste0('functions/rule.process.R'), encoding = 'UTF-8')
source(paste0('functions/CRL_functions.R'), encoding = 'UTF-8')

path <- 'simulation/study one/'
data.path <- paste0(path,'data/')

# define data settings
Ns <- c(2500,5000,10000)
eff.sizes <- seq(5,30,5)
num.true.grps <- c(1,3,5)

# number of re-sampled data for each setting
num.run <- 100 
num.tree <- 50

# change the variable name of your outcome of interest to 'outcome' and your treatment to 'treat'
# before running the model, ensure your training data and test data are of same form.
for (N in Ns){
  for (eff.size in eff.sizes) {
    for (num.true.grp in num.true.grps) {  
      # define parameters and result container
      if(N==2500){
        min.node.size <- case_when(
          num.true.grp == 1  ~ 70
          , num.true.grp == 3 ~ 50
          , num.true.grp == 5 ~ 80
        )
      }
      if(N==5000){
        min.node.size <- case_when(
          num.true.grp == 1  ~ 160
          , num.true.grp == 3 ~ 80
          , num.true.grp == 5 ~ 150
        )
      }
      if(N==10000){
        min.node.size <- case_when(
          num.true.grp == 1  ~ 300
          , num.true.grp == 3 ~ 150
          , num.true.grp == 5 ~ 80
        )
      }
      
      # load data 
      datafile <- paste0(N,'obs;',eff.size,'tau;',num.true.grp,'grps')
      paste0('This is for data: ',datafile) %>% print()
      load(paste0(data.path,datafile,'.RData'))
      
      # create model result container
      cate.train <- list() # all estimates on training samples
      cate.test <- list() # all estimates on test samples
      rule.bi.trs <- list() # all 0/1 data representing the satisfaction of rules of training samples
      rule.bi.tes <- list() # all 0/1 data representing the satisfaction of rules of test samples
      selection_models <- list() # all D-learning models learnt
      Rules <- list() # all rules generated
      rules <- list() # all rules selected
      Subgrps <- list() # all rules after analysis
      FDRs <- list() # FDR results in rule significance analysis
      
      # run model
      for(i in 1:num.run){
        #i <- 1
        print(i)
        train <- data.train[[i]]
        test <- data.test[[i]]
        # you have to provide your (estimated) propensity scores based on your case
        # in our case, samples are purely randomized.
        ps.train <- rep(0.5,nrow(train))
        ps.test <- rep(0.5,nrow(test))
        
        result_discovery <- CRL_rule_discovery(train = train,test = test,
                                               name.covars = paste0('x',1:9),
                                               min.node.size = 40,num.tree = 50,
                                               seed = 2025)
        # save discovery results 
        rule.bi.trs[[i]] <- result_discovery$rule.bi.tr
        rule.bi.tes[[i]] <- result_discovery$rule.bi.te
        Rules[[i]] <- result_discovery$Rules
        
        result_selection <- CRL_rule_selection(train = train,
                                               rule.bi.tr = result_discovery$rule.bi.tr,
                                               Rules = result_discovery$Rules,
                                               ps.train = ps.train        
                                               )
        # save selection results 
        selection_models[[i]] <- result_selection$model
        rules[[i]] <- result_selection$rules_selected
        
        
        #---- estimate HTE for samples in training and test data with the learnt D-learning model ----
        cate.train[[i]] <- predict(result_selection$model,newx = as.matrix(result_discovery$rule.bi.tr[-1]))
        cate.test[[i]] <- predict(result_selection$model, newx = as.matrix(result_discovery$rule.bi.te[-1]))
        
        result_analysis <- subgrp.analysis(result_selection$rules_selected,
                                           compare = cate.train[[i]],
                                           rule.bi = result_discovery$rule.bi.tr,
                                           df = train,
                                           alpha = 0.05,
                                           FDR = 0.05
        )
        
        # save analysis results
        Subgrps[[i]] <-  result_analysis$fin.subgrp
        FDRs[[i]] <- result_analysis$FDR
      }
    }
    
    # save all results
    resultfile <- paste0(datafile,'_result')
    save(Rules,rules,selection_models,Subgrps,cate.train,cate.test,rule.bi.trs,rule.bi.tes,FDRs,
         file = paste0(path,'CRL/model.result/',resultfile,'.RData'))
  }
}













