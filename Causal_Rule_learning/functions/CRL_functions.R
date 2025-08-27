#' Causal Rule Learning - Rule Discovery Step
#'
#' This function implements the rule discovery step in the Causal Rule Learning (CRL) framework.
#' It first trains a causal forest on the training dataset to estimate heterogeneous treatment effects,
#' then decomposes each tree into candidate rules, removes duplicates, and converts the rules into
#' binary compliance indicators for both training and test sets. Rules with extreme support (too rare
#' or too common) are further removed to improve interpretability and stability.
#'
#' @param train A data.frame containing the training data. Must include:
#'   - outcome (numeric): observed outcome variable
#'   - treat (binary: 0/1): treatment assignment
#'   - covariates: columns specified by `name.covars`
#' @param test A data.frame containing the test data, with the same variable structure as `train`.
#' @param name.covars A character vector specifying the names of covariates used in the causal forest.
#' @param min.node.size Minimum node size for each causal tree (default = 40).
#' @param num.tree Number of trees to grow in the causal forest (default = 50).
#' @param seed Random seed for reproducibility (default = 2025).
#'
#' @return A list with three elements:
#'   - `rule.bi.tr`: Training data with binary indicators of rule compliance,
#'     where each compliance indicator is recoded to its subgroup CATE estimate.
#'   - `rule.bi.te`: Test data with binary indicators of rule compliance,
#'     aligned with the training rules.
#'   - `Rules`: A data.frame of discovered rules and their estimated CATEs.
#'
#' @details
#' The workflow includes:
#' 1. Train a causal forest (`grf::causal_forest`) on the training set.  
#' 2. Decompose each tree into a set of human-readable rules (`grf_tree2rule`).  
#' 3. Remove duplicated or redundant rules via `rule.remove`.  
#' 4. Convert rules into binary indicators (`raw2rule.binary`) for both train and test sets.  
#' 5. Recode compliance indicators: `1` is replaced with the rule’s estimated CATE.  
#' 6. Filter out rules with extreme support (obey proportion < 0.1 or > 0.9).  
#'
#' This step corresponds to the "rule discovery" stage in CRL, producing candidate subgroup rules
#' that are interpretable and clinically meaningful for downstream selection and analysis.
#'
#' @examples
#' # Example (using synthetic data)
#' set.seed(123)
#' train <- data.frame(
#'   outcome = rnorm(200),
#'   treat = sample(0:1, 200, replace = TRUE),
#'   x1 = rnorm(200),
#'   x2 = rnorm(200)
#' )
#' test <- data.frame(
#'   outcome = rnorm(100),
#'   treat = sample(0:1, 100, replace = TRUE),
#'   x1 = rnorm(100),
#'   x2 = rnorm(100)
#' )
#' res <- CRL_rule_discovery(train, test, name.covars = c("x1","x2"), num.tree = 10)
#' str(res)
#'
#' @export
CRL_rule_discovery <- function(train = train, test = test,
                               name.covars,
                               min.node.size = 40, num.tree = 50,
                               seed = 2025){
  
  # Train causal forest
  tau.forest <- grf::causal_forest(
    X = train[name.covars],
    Y = train$outcome, 
    W = train$treat, 
    min.node.size = min.node.size, 
    num.trees = num.tree, 
    seed = seed
  )
  
  # Extract rules from each tree
  list.rules <- list()
  for (j in 1:num.tree){
    t <- grf::get_tree(tau.forest, j)
    list.rules[[j]] <- grf_tree2rule(
      tree = t,
      var.names = name.covars,
      data = train
    )
  }
  
  if(length(list.rules) == 0){
    paste0('Empty rule list from the dataset') %>% print()
    return(NULL)
  }
  
  # Decompose trees into rules and remove duplicates
  Rules <- rule.remove(list.rules, remove.lens = c(), remove.duplicated = TRUE)
  
  # Convert rules into binary indicators
  rule.bi.tr <- raw2rule.binary(data = train, Rules = Rules$rule, keep.vars = 'None')
  rule.bi.te <- raw2rule.binary(data = test, Rules = Rules$rule, keep.vars = 'None')
  
  # Compute rule support (proportion of training samples obeying each rule)
  obey.prop <- c()
  for (j in 1:(ncol(rule.bi.tr)-1)){
    var <- paste0('obey_rule', j)
    obey.prop[j] <- round(sum(rule.bi.tr[var] == 1) / nrow(rule.bi.tr), digits = 2) 
    
    # Recode compliance: 1 -> corresponding rule’s CATE estimate
    rule.bi.tr[var] <- ifelse(rule.bi.tr[var] == 1, Rules$cate[j], 0)
    rule.bi.te[var] <- ifelse(rule.bi.te[var] == 1, Rules$cate[j], 0)
  }
  
  summary(obey.prop)
  
  # Filter rules with extreme support
  if(length(which(obey.prop < 0.1 | obey.prop > 0.9)) != 0){
    rule.bi.tr <- rule.bi.tr[-which(obey.prop < 0.1 | obey.prop > 0.9)]
    rule.bi.te <- rule.bi.te[colnames(rule.bi.tr)]
  }
  
  # Return results
  result <- list(
    rule.bi.tr = rule.bi.tr,
    rule.bi.te = rule.bi.te,
    Rules = Rules
  )
  
  return(result) 
}


#' Causal Rule Learning - Rule Selection Step
#'
#' This function implements the rule selection step in the Causal Rule Learning (CRL) framework.
#' It applies D-learning with LASSO regularization to select the most predictive rules 
#' from the discovered candidate rules. Cross-validation is used to tune the regularization 
#' parameter \eqn{\lambda}, ensuring robust and interpretable subgroup selection.
#'
#' @param train A data.frame containing the training data with columns:
#'   - outcome (numeric): observed outcome variable  
#'   - treat (binary: 0/1): treatment assignment  
#' @param rule.bi.tr A data.frame containing binary indicators (or CATE-coded indicators)
#'   of rule compliance for the training set, as returned by \code{CRL_rule_discovery}.
#' @param Rules A data.frame of discovered rules and their estimated CATE values, 
#'   as returned by \code{CRL_rule_discovery}.
#' @param ps.train A numeric vector of estimated propensity scores for the training set.  
#' 
#' @return A list with the following elements:
#'   - `model`: The fitted \code{glmnet} model object.  
#'   - `rules_selected`: A data.frame of rules that were selected (non-zero coefficients).  
#'   - `rule.ind`: Integer indices of the selected rules.  
#'
#' @details
#' The selection process includes:
#' 1. Constructing the D-learning outcome: \eqn{Y* = 2 * Y * W},  
#'    where \eqn{Y} is the outcome and \eqn{W} is the treatment indicator.  
#' 2. Using rule-based features \eqn{X} from `rule.bi.tr` as predictors.  
#' 3. Weighting each sample by the inverse propensity score \eqn{1 / e(X)}.  
#' 4. Running LASSO-penalized regression via \code{glmnet}, with \eqn{\lambda} tuned by cross-validation.  
#' 5. Selecting rules with non-zero coefficients, representing informative subgroups.  
#'
#' This step corresponds to the "rule selection" stage in CRL, filtering 
#' interpretable rules that best explain treatment effect heterogeneity.
#'
#' @examples
#' # Example (synthetic)
#' set.seed(123)
#' train <- data.frame(
#'   outcome = rnorm(200),
#'   treat = sample(0:1, 200, replace = TRUE)
#' )
#' rule.bi.tr <- data.frame(id = 1:200, obey_rule1 = rbinom(200, 1, 0.3), obey_rule2 = rbinom(200, 1, 0.5))
#' Rules <- data.frame(rule = c("x1 > 0", "x2 <= 5"), cate = c(1.2, -0.8))
#' ps.train <- runif(200, 0.2, 0.8)
#' ps.test <- runif(200, 0.2, 0.8)
#'
#' res <- CRL_rule_selection(train, rule.bi.tr, Rules, ps.train, ps.test)
#' str(res$rules_selected)
#'
#' @export
CRL_rule_selection <- function(train,
                               rule.bi.tr,
                               Rules,
                               ps.train
){
  
  # Prepare data for D-learning
  y.new <- 2 * train$outcome * train$treat
  X <- as.matrix(rule.bi.tr[-1])
  weights <- 1 / ps.train 
  
  # Cross-validation to tune lambda
  cv.rule.model <- tryCatch({
    glmnet::cv.glmnet(
      x = X,
      y = y.new,
      weights = weights,
      type.measure = 'mse',
      nfolds = 10,
      trace.it = 1
    )
  },
  error = function(e){
    message('Error in cross-validation');
    return(NULL)
  })
  
  # Fit model with best lambda
  d.rule.model <- glmnet::glmnet(
    x = X,
    y = y.new,  
    family = 'gaussian',
    weights = weights,
    alpha = 1,  # LASSO penalty
    lambda = cv.rule.model$lambda.min,
    standardize = TRUE,
    intercept = FALSE,
    trace.it = TRUE
  )
  
  # Identify non-zero (selected) rules
  ind.which <- which(as.vector(d.rule.model$beta) != 0)
  d.rule.codes <- colnames(X)[ind.which] %>%
    gsub(pattern = 'obey_rule', replacement = '') %>%
    as.integer()
  
  rules <- Rules[d.rule.codes,]
  
  # Return results
  result <- list(
    model = d.rule.model,
    rules_selected = rules,
    rule.ind = d.rule.codes
  )
  
  return(result) 
}



#' Perform Subgroup Analysis with False Discovery Rate (FDR) Control available
#'
#' This function analyzes subgroups (rules) from the rule selection results, performing
#' significance testing on both complete rules and their decomposed versions while
#' controlling the False Discovery Rate. It helps identify robust subgroups that
#' show significant treatment effect differences.
#'
#' @param subgrp A data.frame containing subgroups (rules) from a causal inference method.
#' Should contain a 'rule' column with rule definitions and have row names
#' that correspond to rule indices.
#' @param compare A numeric vector of CATE (Conditional Average Treatment Effect)
#' estimates for each sample in the dataset.
#' @param rule.bi A data.frame containing binary rule endorsement indicators for
#' each sample (typically output from raw2rule.binary).
#' @param df The original raw data used to evaluate the rules.
#' @param alpha Significance level for individual rule testing. Default is 0.05.
#' @param FDR False Discovery Rate threshold for multiple testing correction.
#' Set to FALSE to disable FDR correction. Default is 0.05.
#'
#' @return A list with two components:
#' \describe{
#' \item{fin.subgrp}{A refined data.frame of subgroups containing only rules
#' that passed both complete rule significance testing and
#' decomposition analysis.}
#' \item{FDR}{If FDR != FALSE, a data.frame with FDR testing results including
#' original p-values, FDR-adjusted p-values, and inclusion decision
#' for each rule.}
#' }
#' If only one rule exists in subgrp, returns the original rule without analysis.
#'
#' @export
#'
#' @examples
#' # Example usage (assuming appropriate data structures):
#' # results <- subgrp.analysis(
#' # subgrp = rules_df,
#' # compare = cate_estimates,
#' # rule.bi = rule_endorsement_matrix,
#' # df = original_data,
#' # alpha = 0.05,
#' # FDR = 0.05
#' # )
#'
#' # View refined subgroups
#' # print(results$fin.subgrp)
#'
#' # View FDR results (if applicable)
#' # print(results$FDR)

subgrp.analysis <- function(subgrp, compare, rule.bi, df, alpha = 0.05, FDR = 0.05){
  res <- list()
  if(nrow(subgrp) <= 1){fin.subgrp <- subgrp} # no anlysis if only 1 group exists
  else{ 
    #- significance analysis (complete rule)-
    rule.inds <- row.names(subgrp)
    endorse.cols <- rule.bi[paste0('obey_rule',rule.inds)]
    
    pvals.comp.ori <- pvals <- lapply(endorse.cols, FUN = rule.analysis,compare = compare) %>% unlist()
    
    if(FDR!=FALSE){ 
      #use FDR control for rule significance analysis, instead of the original p-value
      res[['FDR']] <- subgrp
      res[['FDR']]$p_ori <- pvals.comp.ori
      res[['FDR']]$p_fdr <- pvals <- p.adjust(pvals, method = "BH")
      res[['FDR']]$res <- ifelse(res[['FDR']]$p_fdr<=FDR,'in','out')
    }
    rule.inds <- rule.inds[pvals < alpha]
    #- decomposition analysis -
    decomp.pvals <- list() #save the p-values of the decomposed rules
    decomp.binary <- data.frame(id = 1:nrow(df))
    for (i in rule.inds){
      #print(i)
      rule <- subgrp[i,'rule'] # get a specific rule
      # decompose the current rule
      rule <- strsplit(x = rule, split = ' & ', fixed = T)[[1]]
      num.cond <- length(rule) 
      if(num.cond > 1){ # need to be decomposed 
        temp <- c() 
        for (p in 1:num.cond) { # p refers to the p-th factor to be removed    
          # remove the p-th factor in the rule
          # p <- 2
          col.name <- paste0('obey',i,'_remove_',p)
          new.rule <- paste0(rule[c(1:num.cond)[-p]],collapse = ' & ') 
          # determine whether the revised rule is endorsed
          command <- paste0('as.factor(ifelse(',new.rule,', 1 ,0))') 
          decomp.binary[col.name] <- with(df,{parse(text = command) %>% eval()})
          temp[p] <- rule.analysis(compare = compare, endorse.col = decomp.binary[col.name])
        }
        decomp.pvals[[i]] <- temp}
      else{decomp.pvals[[i]] <- 10000} # a big number to ensure length-one rule pass the test
    }
    #- delete rules according to the p-values - 
    # all the revised rules should have larger or equal p-values than that of the complete rule
    flags <- lapply(rule.inds, function(i){
      ifelse(min(decomp.pvals[[i]]) >= pvals.comp.ori[paste0('obey_rule',i)],T,F)}) %>% unlist()
    res[['fin.subgrp']] <- subgrp[flags,] 
  }
  return(res)
}



