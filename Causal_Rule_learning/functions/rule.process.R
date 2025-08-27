

#' Decompose a Causal Tree into Rules (Subgroups) with CATE Estimates
#'
#' This function takes a causal tree (e.g., from a grf object) and decomposes it into a set of rules.
#' Each rule defines a subgroup based on the splitting criteria along the path from the root to a leaf node.
#' For each subgroup (leaf node), the function calculates the corresponding Conditional Average Treatment Effect (CATE).
#'
#' @param tree A tree object (list structure), typically extracted from a grf forest using get_tree().
#' @param var.names A character vector containing the names of the features/variables used in the tree.
#' @param data A data.frame containing the original data used to grow the tree. Must include columns named 'treat' (treatment indicator) and 'outcome'.
#'
#' @return A data.frame with rows corresponding to each leaf node (subgroup) and the following columns:
#' \describe{
#' \item{node}{Integer. The node ID of the leaf in the original tree structure.}
#' \item{rule}{Character. The rule defining the subgroup, constructed from the split conditions on the path from the root to this leaf.}
#' \item{cate.hat}{Numeric. The estimated Conditional Average Treatment Effect (CATE) for the subgroup defined by the rule.}
#' }
#' Returns NULL if the input tree consists of only a root node (i.e., no splits).
#' @export
#'
#' @examples
#' # Assuming tau.forest is a trained causal_forest object from grf
#' # and X is the matrix of covariates used to train it
#' # tree <- get_tree(tau.forest, 1)
#' # rule_df <- grf_tree2rule(tree = tree,
#' # var.names = colnames(X),
#' # data = data.frame(X, treat = treatment_vector, outcome = outcome_vector))

grf_tree2rule <- function(tree = tree, var.names = var.names, data = data){
  # Reorganize tree structure and create a simplified tree structure 
  nodes <- tree$nodes
  if(length(nodes)==1) return(NULL)
  else{
    df.nodes <- data.frame(node = 1:length(nodes), 
                           father = c(NA,1,1,rep(NA,length(nodes)-3)),
                           left.child = 2,
                           right.child = 3)
    for (i in 2:length(nodes)){ 
      #traverse all nodes
      if(nodes[[i]]$is_leaf){ #define the child to be NA for leaf nodes
        df.nodes[i,'left.child'] <- NA
        df.nodes[i,'right.child'] <- NA
      }
      else{ #define the left and right child and their parent for internal nodes
        df.nodes[i,'left.child'] <- nodes[[i]]$left_child
        df.nodes[i,'right.child'] <- nodes[[i]]$right_child
        df.nodes[nodes[[i]]$left_child,'father'] <- df.nodes[nodes[[i]]$right_child,'father'] <- i
      }
    }
    
    #---- extract internal node components ----
    components <- list() 
    
    # extract components from all internal nodes
    for (i in 1:length(nodes)){
      if(!nodes[[i]]$is_leaf){ 
        split.var <- var.names[nodes[[i]]$split_variable]
        split.val <- round(nodes[[i]]$split_value, digits = 2)
        # extract the left and right components 
        com.left <- paste0(split.var,' <= ',split.val)
        com.right <- paste0(split.var,' > ',split.val) 
        # save components
        components[[as.character(i)]]<- c(com.left,com.right)}
    }  
    
    # obtain the index for internal nodes and leaf nodes
    inter.nodes <- names(components) %>% as.integer()
    leaf.nodes <- setdiff(1:length(nodes),inter.nodes)
    
    #---- get the path to each leaf node of the tree (bottom-up)----
    for (leaf in leaf.nodes) {
      # print(leaf)
      rule <- c()  #initialize a tree
      cur.node <- leaf
      while (cur.node>1) {#stop when reach the root
        father.node <- df.nodes[cur.node,'father']
        if(df.nodes[father.node,'left.child']== cur.node){
          comp <- components[[as.character(father.node)]][1]}else{
            comp <- components[[as.character(father.node)]][2]}
        
        rule <- c(comp,rule)
        cur.node <- father.node
      }
      # save the current rule 
      df.nodes[leaf,'rule'] <- paste(rule,collapse = ' & ')
      #print(df.nodes[leaf,'rule'])
    }
    
    #preserve leaf nodes only
    df.nodes <- df.nodes[leaf.nodes,c('node','rule')]
    row.names(df.nodes) <- 1:length(leaf.nodes)
    
    # simplify and standardize the rules
    df.nodes$rule <- lapply(df.nodes$rule, 
                            function(rule){rule.simp(rule)}) %>% unlist
    
    # for each leaf node, calculate its correspoonding cate estimates
    for(leaf in df.nodes$node){
      #print(leaf)
      samples <-  nodes[[leaf]]$samples
      treat.samples <- samples[which(data[samples,'treat']==1)] 
      control.samples <- setdiff(1:nrow(data),treat.samples)
      df.nodes[df.nodes$node == leaf,'cate.hat'] <- mean(data[treat.samples,'outcome'])- mean(data[control.samples,'outcome'])
    }
    return(df.nodes) 
  }
}

#' Round Threshold Values in a Rule String
#'
#' This function takes a rule string (e.g., "age <= 25.6 & income > 40000.78"),
#' extracts the numeric threshold values, rounds them to the specified number of digits,
#' and reconstructs the rule string with the rounded values.
#'
#' @param rule A character string representing a rule, composed of conditions joined by ' & '.
#' Each condition should be in the format "variable operator value" (e.g., "age <= 25.6").
#' @param digit An integer specifying the number of decimal places to round the numeric values to. Default is 0 (round to whole numbers).
#'
#' @return A character string of the reconstructed rule with rounded threshold values.
#' For example, input "age <= 25.6 & income > 40000.78" with digit = 0 returns "age <= 26 & income > 40001".
#'
#' @export
#'
#' @examples
#' rule <- "age <= 25.6 & income > 40000.78"
#' rule.thre.rd(rule, digit = 0) # Returns "age <= 26 & income > 40001"
#' rule.thre.rd(rule, digit = 1) # Returns "age <= 25.6 & income > 40000.8"


rule.thre.rd <- function(rule, digit = 0){
  # Split the rule string into individual components/conditions
  comps <- strsplit(rule, split = ' & ')[[1]]
  # Initialize a dataframe to store parsed components
  df.rule <- data.frame(var = NA, sign = NA, value = NA)
  # Parse each component into variable, operator, and value
  for (i in 1:length(comps)){
    df.rule[i,] <- stringr::str_match(comps[i], pattern="([[:alnum:]]+)([<>=]+)([-\\d.]+)")[2:4]
  }
  df.rule$value <- as.numeric(df.rule$value)
  df.rule$value <- round(df.rule$value, digits = 0)
  # recombine the components
  for (i in 1:length(comps)) {
    comps[i] <- paste0(df.rule[i,], collapse = '')
  }
  final.rule <- paste(comps, collapse  = ' & ')
  return(final.rule)
}
#' Adjust Threshold Value for a Specific Variable in a Rule String
#'
#' This function searches for conditions involving a specified variable within a rule string.
#' If the variable's threshold value is within an acceptable range of deviation from a desired true threshold,
#' it replaces the value with the true threshold. Otherwise, the original value is kept.
#'
#' @param rule A character string representing a rule, composed of conditions joined by ' & '.
#' Each condition should be in the format "variable operator value" (e.g., "x3 <= 4.91").
#' @param var.name A character string specifying the name of the variable to check and potentially adjust. Default is 'x3'.
#' @param var.thre A numeric value representing the desired "true" threshold to use for adjustment.
#' @param acc.ran A numeric value specifying the acceptable range of deviation from var.thre.
#' If the absolute difference between the existing value and var.thre is <= acc.ran,
#' the value will be replaced with var.thre. Default is 0.1.
#'
#' @return A character string of the reconstructed rule. The threshold for the specified variable is adjusted
#' if it falls within the acceptable range of the desired threshold; otherwise, the original rule is returned.
#' If the specified variable does not appear in the rule, the original rule is returned unchanged.
#'
#' @export
#'
#' @examples
#' rule <- "x1 > 2.5 & x3 <= 4.91 & x2 != 0"
#' # Adjust 'x3' if its value is within 0.1 of 5.0
#' rule.thre.adjust(rule, var.name = 'x3', var.thre = 5, acc.ran = 0.1)
#' # Returns: "x1 > 2.5 & x3 <= 5 & x2 != 0"
#'
#' # Variable 'x4' not in rule, so original rule is returned
#' rule.thre.adjust(rule, var.name = 'x4', var.thre = 10, acc.ran = 1)
#' # Returns: "x1 > 2.5 & x3 <= 4.91 & x2 != 0"

rule.thre.adjust <- function(rule, var.name = 'x3', 
                             var.thre = 5, acc.ran = 0.1){
  # rule <- "x3 <= 4.91" 
  # decompose the rule
  df.rule <- data.frame(var = NA, sign = NA, value = NA)
  comps <- strsplit(rule, split = ' & ')[[1]]
  if(length(comps)==1){
    df.rule[i,] <- strsplit(comps, split = ' ')[[1]]
  }else{
    for (i in 1:length(comps)){
      df.rule[i,] <- stringr::str_match(comps[i], pattern="([[:alnum:]]+)([<>=]+)([-\\d.]+)")[2:4]
    }
  }
  df.rule$value <- as.numeric(df.rule$value)
  inds <- which(df.rule$var == var.name)
  if(length(inds)== 0) return(rule)
  else{
    # if the cut-off value is within the given range,
    # reassign it to be the desired threshold,
    # otherwise, keep its original value
    df.rule[inds,'value'] <- lapply(df.rule[inds,'value'],
                                    function(value){
                                      ifelse(abs(value-var.thre)<= acc.ran, 
                                             var.thre,value)})  
  }
  # recombine the components
  for (i in 1:length(comps)) {
    comps[i] <- paste0(df.rule[i,], collapse = '')
  }
  final.rule <- paste(comps, collapse  = ' & ')
  return(final.rule)
}

#' Round Cut-off Points in a Rule String to Integer or Decimal Places
#'
#' This function processes a rule string and rounds the numeric threshold values
#' based on the variable names. Specific variables are rounded to integers (0 decimal places),
#' while all other variables are rounded to 2 decimal places.
#'
#' @param rule A character string representing a rule, composed of conditions joined by ' & '.
#' Each condition should be in the format "variable operator value" (e.g., "sysbp <= 120.5").
#' @param from.me A character string indicating the source or format of the rule.
#' Currently supports 'CRL' and 'CRE' formats (standard pattern matching)
#' or 'PRIM' format (character splitting). Default is 'CRL'.
#'
#' @return A character string of the reconstructed rule with threshold values rounded appropriately.
#' Variables in the predefined list are rounded to integers, others to 2 decimal places.
#'
#' @export
#'
#' @examples
#' rule <- "sysbp <= 120.51 & age > 45.67 & bmi < 25.123"
#' rule.cp2integer(rule, from.me = 'CRL')
#' # Returns: "sysbp <= 121 & age > 46 & bmi < 25.12"
#' 
rule.cp2integer <-  function(rule, from.me = 'CRL'){
  comps <- strsplit(rule, split = ' & ')[[1]]
  df.rule <- data.frame(var = NA, sign = NA, value = NA)
  for (i in 1:length(comps)){
    if(from.me %in% c('CRL','CRE')){
      df.rule[i,] <- stringr::str_match(comps[i], pattern="([[a-zA-Z0-9.]]*)([<>=]+)([-\\d.]+)")[2:4]
    }
    if(from.me == 'PRIM'){
      df.rule[i,] <-  strsplit(comps[i], split = '')[[1]]
    }
  }
  # round or keep two digits according to specific variables
  which.ind <- df.rule$var %in% c('sysbp','age','othdhis','chdhis','othchd','atrial',
                                  'nyha','mesh','phth','murmur','electb','xece',
                                  'diabp','sex','surgesitua','oroom','occnum','periop','employ')
  
  df.rule$value[which.ind] <- round(as.numeric(df.rule$value[which.ind]), digits = 0)
  df.rule$value[!which.ind] <- round(as.numeric(df.rule$value[!which.ind]), digits = 2)
  
  # recombine the components
  for (i in 1:length(comps)) {
    comps[i] <- paste0(df.rule[i,], collapse = '')
  }
  final.rule <- paste(comps, collapse  = ' & ')
  return(final.rule)
}

#' Simplify and Standardize a Rule by Merging Redundant Conditions
#'
#' This function takes a rule string containing multiple conditions joined by ' & ',
#' identifies conditions on the same variable with the same operator, and merges them
#' by taking the most restrictive value (max for '>' operators, min for '<' operators).
#'
#' @param rule A character string representing a rule, composed of conditions joined by ' & '.
#' Each condition should be in the format "variable operator value" (e.g., "age <= 30 & age > 20").
#'
#' @return A character string of the simplified rule with redundant conditions merged.
#' If multiple conditions exist for the same variable-operator combination,
#' they are replaced with a single condition using the most restrictive value.
#' For a single-condition rule, the original rule is returned unchanged.
#'
#' @export
#'
#' @examples
#' rule <- "age > 20 & age > 25 & income <= 50000 & income <= 60000"
#' rule.simp(rule)
#' # Returns: "age > 25 & income <= 50000"
#'
#' rule2 <- "age <= 30"
#' rule.simp(rule2) # Returns: "age <= 30"
#'
rule.simp <- function(rule){
  rule <- strsplit(rule,split = ' & ')[[1]]
  # for length-one rule, no need to simplify but return 
  if(length(rule)==1){final.rule <- rule} 
  else{ 
    df.rule <- data.frame(var = NA,sign = NA,value = NA)
    sta <- strsplit(rule, split = ' ') 
    for (i in 1:length(sta)) {df.rule[i,] <-  sta[[i]]}
    df.rule$var_sign <- paste0(df.rule$var,df.rule$sign)
    # print(df.rule)
    new.rule <- c()
    for(each in unique(df.rule$var_sign)){
      values <- df.rule[which(df.rule$var_sign == each),'value'] 
      new.comp <- paste0(each,
                         ifelse(strsplit(each,split = '')[[1]][nchar(each)]=='>',
                                max(values),min(values)))
      
      new.rule <- c(new.rule,new.comp)
    }
    final.rule <- paste(new.rule, collapse  = ' & ')
  }
  return(final.rule)
}


#' Format CRL Rule Expressions
#'
#' This function reformats a rule string used in Causal Rule Learning (CRL).
#' It replaces inequality-based conditions (`>0`, `<=0`) with equivalent
#' binary indicator expressions (`==1`, `==0`), making the rules easier to 
#' interpret in a binary feature space.
#'
#' @param rule A character string representing a rule, e.g. 
#'   `"x1<=0 & x2>0 & x3<=5"`.
#'
#' @return A character string with reformatted rule expressions.
#'
#' @examples
#' # Example of reformatting a rule
#' format_CRL_rule("x1<=0 & x2>0 & x3<=5")
#' # Output: "x1==0 & x2==1 & x3<=5"
#'
#' @export
format_CRL_rule <- function(rule){
  # Replace >0 with ==1, and <=0 with ==0
  rule <- gsub('>0','==1',rule)
  rule <- gsub('<=0','==0', rule)
  return(rule)
}


#' Format Causal Tree (CT) Rules for Comparison with CRL
#'
#' This function reformats rules generated by a Causal Tree (CT) model
#' into a unified format so that they can be compared directly with 
#' Causal Rule Learning (CRL) rules.
#'
#' @param df.rules A data.frame containing CT rules, where:
#'   - `outcome` column stores the estimated treatment effect (CATE).
#'   - Other columns store rule conditions for each variable.
#'
#' @return A data.frame with:
#'   - `rule`: character strings representing formatted rules  
#'   - `cate`: numeric estimated treatment effects
#'
#' @details
#' The formatting steps include:
#' - Replacing `"is"` with `"=="`  
#' - Trimming whitespace and removing empty elements  
#' - Sorting conditions by variable name  
#' - Adjusting inequality operators (`>=` → `>`, `<` → `<=`)  
#' - Collapsing conditions into a single string with `" & "`  
#'
#' If no rules are found, the function returns a data frame with `NA` in `rule`.
#'
#' @examples
#' # Example input resembling CT rule output
#' df.rules <- data.frame(
#'   outcome = c(1.2, -0.8),
#'   x1 = c("x1 is 1", "x1 is 0"),
#'   x2 = c("x2 < 5", "")
#' )
#' format_CT_rule(df.rules)
#'
#' @export
#'
format_CT_rule <- function(df.rules){
  df <- data.frame(rule = NA, cate = df.rules$outcome)
  
  if(!is.na(df$rule[1])){
    rules <- df.rules[3:ncol(df.rules)]
    for (p in 1:nrow(rules)){
      r <- rules[p,]
      # Replace "is" with "=="
      r <- gsub('is','==',r)
      # Remove extra spaces
      r <- trimws(r) 
      r <- r[which(r!="")] 
      # Collapse into a single rule string
      r <- paste0(r, collapse = '')
      # Sort conditions by variable name
      r <- strsplit(r, split = '&')[[1]]
      r <- paste(sort(r), collapse = ' & ')
      # Adjust inequality formats
      r <- gsub('>=','>',r)
      r <- gsub('<','<=',r)
      # Store formatted rule
      df[p,'rule'] <- r 
    }
  }
  
  # If no rules are extracted, return NA
  return(df)
}

#' Remove Duplicated Rules and Aggregate Estimates
#'
#' This function removes duplicated rules from a data frame and replaces
#' other columns (e.g., CATE estimates, coefficients) with aggregated values.
#' For duplicated rules, the `cate.hat` values are averaged, and if a `coef`
#' column exists, the coefficients are updated proportionally to the averaged
#' treatment effect.
#'
#' @param temp A data.frame containing:
#'   - `rule` (character): rule string  
#'   - `cate.hat` (numeric): estimated CATE  
#'   - `coef` (numeric, optional): coefficient/weight for the rule
#'
#' @return A data.frame with:
#'   - Unique rules (`rule`)
#'   - Aggregated `cate.hat` values (mean across duplicates)
#'   - Updated `coef` values if present  
#'   Row names are reset to sequential indices.
#'
#' @details
#' - For each duplicated rule, `cate.hat` values are averaged across duplicates.  
#' - If a `coef` column exists, the coefficient is recalculated as a weighted sum 
#'   of CATE estimates divided by the averaged CATE.  
#' - All duplicated rows except the first occurrence are removed.  
#'
#' @examples
#' # Example with duplicated rules
#' temp <- data.frame(
#'   rule = c("x1 > 0 & x2 < 5", "x1 > 0 & x2 < 5", "x3 > 10"),
#'   cate.hat = c(1.2, 1.4, -0.8),
#'   coef = c(0.5, 0.7, 0.3)
#' )
#' rm.duplicates(temp)
#'
#' @export
rm.duplicates <- function(temp){
  duplicates <- temp[duplicated(temp$rule),]
  
  # For each duplicated rule, average cate.hat and update coef if present
  if(nrow(duplicates)!=0){
    for(i in 1:nrow(duplicates)){
      dup <- duplicates$rule[i][[1]]
      obs <- which(temp$rule == dup)
      
      # Always average the CATE estimate
      temp[obs[1],'cate.hat'] <- mean(temp[obs, 'cate.hat'], na.rm = TRUE)
      
      # If coef column exists, recompute weight after cate.hat averaging
      if ('coef' %in% colnames(temp)){ 
        weighted.sum.cate <- temp[obs,'cate.hat'] %*% temp[obs,'coef']
        temp[obs[1],'coef'] <- weighted.sum.cate / temp[obs[1],'cate.hat']
      }
    }
    temp <- temp[!duplicated(temp$rule),]
  }
  
  row.names(temp) <- 1:nrow(temp) 
  return(temp)
}
#' Remove Duplicated or Specified-Length Rules
#'
#' This function processes a list of rules by removing those with specified lengths
#' and eliminating duplicated rules. It also recalculates subgroup CATE (Conditional
#' Average Treatment Effect) estimates for the remaining unique rules.
#'
#' @param rulelist A list of data frames, each containing at least:
#'   - `rule` (character): rule string
#'   - `cate.hat` (numeric): estimated CATE
#' @param remove.lens A numeric vector specifying rule lengths (number of conditions 
#'   joined by "&") to remove. Default is an empty vector.
#' @param remove.duplicated Logical, whether to remove duplicated rules. 
#'   Default = TRUE.
#'
#' @return A data frame with columns:
#'   - `rule`: unique rule strings (character)
#'   - `cate.hat`: corresponding CATE estimates (numeric)  
#'   Returns `NULL` if no rules remain after filtering.
#'
#' @details
#' - Rules are reorganized so that their components are sorted consistently 
#'   (e.g., x1, x2, ...).
#' - Duplicate removal relies on an auxiliary function `rm.duplicates()`.
#' - A message is printed if duplicate removal is successful.
#'
#' @examples
#' # Example: remove rules of length 2 and eliminate duplicates
#' rulelist <- list(
#'   data.frame(rule = "x1 > 0 & x2 < 5", cate.hat = 1.2),
#'   data.frame(rule = "x2 < 5 & x1 > 0", cate.hat = 1.3),
#'   data.frame(rule = "x3 > 10", cate.hat = -0.8)
#' )
#' rule.remove(rulelist, remove.lens = 2, remove.duplicated = TRUE)
#'
#' @export
rule.remove <- function(rulelist, remove.lens = c(), remove.duplicated = T){ 
  # Extract rules from list
  rule <- lapply(rulelist,function(df){if(!is.null(df)) df$rule}) %>% unlist %>% as.character()
  cate.hat <- lapply(rulelist,function(df){if(!is.null(df)) df$cate.hat}) %>% unlist
  Rules <- data.frame(rule = rule, cate.hat = cate.hat)
  Rules$rule <- as.character(Rules$rule)
  
  # Compute length of all rules
  rule.lens <- lapply(1:nrow(Rules), function(i){
    length(strsplit(Rules[i,'rule'], split = ' & ')[[1]])} ) %>% unlist()
  
  # Remove rules with specified length
  if(length(remove.lens)!=0){
    for (j in remove.lens){
      Rules[rule.lens==j,'rule'] <- NA 
    }
    Rules <- Rules[!is.na(Rules$rule),]
    if(nrow(Rules)==0) return(NULL)
    row.names(Rules) <- 1:nrow(Rules)
  }
  
  # Remove duplicated rules and recalculate CATE estimates
  if(remove.duplicated){
    # Reorganize rule components in ascending variable order
    P <- lapply(Rules$rule, strsplit,split = ' & ' ) %>% unlist(recursive = F)
    temp <- c()
    for (i in 1:length(P)){
      temp[i] <- paste0(sort(P[[i]]),collapse = ' & ') 
    }
    temp = data.frame(rule = temp, cate.hat = Rules$cate.hat,
                      stringsAsFactors = F)
  }
  
  temp <- rm.duplicates(temp)
  
  # Test if there are no duplicates
  if(remove.duplicated & (duplicated(temp$rule) %>% sum() == 0)) 
    print('Rules successfully removed!')
  
  return(temp)
}


#' Convert Raw Data to Rule Compliance Binary Indicators
#'
#' This function evaluates a set of rules on a dataset and creates binary indicators
#' (0/1) for each rule, showing whether each observation complies with the rule.
#' It can also preserve specified original variables in the output.
#'
#' @param data A data.frame containing the raw data to be evaluated against the rules.
#' @param Rules A character vector of rules to be evaluated. Each rule should be a
#' valid R logical expression that can be evaluated in the context of data
#' (e.g., "age > 30 & income < 50000").
#' @param Rules.index An integer vector specifying indices for the rules. Default is
#' sequential numbers from 1 to the number of rules.
#' @param keep.vars A character vector of variable names from the original data to
#' preserve in the output. Use 'None' if no variables should be kept.
#'
#' @return A data.frame with binary compliance indicators for each rule. The output
#' contains:
#' - Preserved original variables (if keep.vars != 'None')
#' - Binary columns named 'obey_ruleX' where X is the rule index, with values:
#' * 1: Observation complies with the rule
#' * 0: Observation does not comply with the rule
#'
#' @export
#'
#' @examples
#' # Create example data
#' set.seed(123)
#' data <- data.frame(
#' age = rnorm(100, 50, 10),
#' income = rnorm(100, 50000, 10000),
#' treatment = sample(0:1, 100, replace = TRUE)
#' )
#'
#' # Define rules
#' rules <- c("age > 45 & income < 60000", "treatment == 1 & age <= 60")
#'
#' # Create binary rule compliance data, keeping original variables
#' rule_binary_data <- raw2rule.binary(data = data,
#' Rules = rules,
#' keep.vars = c('age', 'income', 'treatment'))
#'
#' # View the result
#' head(rule_binary_data)

raw2rule.binary <- function(data = data, Rules, Rules.index = 1:length(Rules),keep.vars){
  if(keep.vars!='None'){
    data.rule <- as.data.frame(data[keep.vars])
    colnames(data.rule) <- keep.vars
    k <- length(keep.vars) #preserve new columns
  }else{
    data.rule <- data.frame(id = 1:nrow(data))
    k <- 1
  }
  num.col <- k+length(Rules)
  data.rule[seq(k+1,num.col)] <- NA 
  colnames(data.rule)[seq(k+1,num.col)] <- paste0('obey_rule',Rules.index) 
  if(ncol(data.rule)!= k+length(Rules)){ print('mistake!')}
  
  #--use parse(), eval() for a series of command
  for (i in 1:length(Rules.index)){
    # print(i)
    txt = paste0('as.factor(ifelse(',Rules[i],', 1 ,0))')
    data.rule[paste0('obey_rule',Rules.index[i])] <-  with(data, 
                                                           expr = {parse(text = txt) %>% eval() })
  }
  data.rule <- as.data.frame(data.rule)
  #data.rule <- as.numeric(data.rule)
  print('binary rule-based data successfully generated!')
  return(data.rule)
}

#' Analyze Rule Compliance Group Differences using Kolmogorov-Smirnov Test
#'
#' This function performs a Kolmogorov-Smirnov (KS) test to determine whether
#' the distribution of a continuous variable differs significantly between
#' rule-compliant and non-compliant groups. This is useful for assessing if
#' a rule effectively identifies subgroups with different characteristics.
#'
#' @param compare A numeric vector of values to be compared between groups
#' (e.g., CATE estimates, positive outcomes, or any continuous measure).
#' @param endorse.col A vector (typically numeric or logical) indicating rule compliance.
#' Non-zero values indicate compliance, zero values indicate non-compliance.
#'
#' @return A numeric value representing the p-value from the two-sample
#' Kolmogorov-Smirnov test. A small p-value (typically < 0.05) suggests
#' that the distributions differ significantly between compliant and
#' non-compliant groups.
#'
#' @export
#'
#' @examples
#' # Example data: CATE estimates and rule compliance
#' set.seed(123)
#' n <- 100
#' cate_estimates <- rnorm(n, 0.5, 0.2) # Simulated CATE estimates
#' rule_compliance <- sample(0:1, n, replace = TRUE) # Rule compliance (0/1)
#'
#' # Test if CATE estimates differ between compliers and non-compliers
#' p_value <- rule.analysis(compare = cate_estimates,
#' endorse.col = rule_compliance)
#' print(p_value)
#'
rule.analysis <- function(compare, endorse.col){
  # endorse.col = endorse.cols['obey_rule1']
  # compare = d.cate
  
  y.grp <- compare[endorse.col!=0]
  y.nongrp <- compare[endorse.col==0]
  pval <- ks.test(y.grp,y.nongrp)$p.value
  
  return(pval)
}

#' Compute Neyman ATE for a Population Defined by Rule Compliance
#'
#' This function calculates the Neyman estimate of the Average Treatment Effect (ATE)
#' for a specific subpopulation defined by their compliance with a rule. The subpopulation
#' can be either those who obey the rule (endorse.col != 0) or those who disobey the rule (endorse.col == 0).
#'
#' @param endorse.col A vector (typically numeric or logical) indicating rule compliance.
#' Non-zero values typically indicate compliance, zero values indicate non-compliance.
#' @param outcome A numeric vector of observed outcomes for the entire population.
#' @param treat A numeric vector of treatment assignments (typically 1 for treatment, 0 for control)
#' for the entire population.
#' @param group A character string specifying which subpopulation to analyze.
#' Options: 'obey' (rule compliers) or 'disobey' (rule non-compliers).
#' Default is 'obey'.
#' @param propensity An optional numeric vector of propensity scores for inverse probability weighting.
#' If NULL (default), no weighting is applied (equivalent to all weights = 1).
#' If provided, should be the same length as outcome and treat.
#'
#' @return A numeric value representing the Neyman estimate of the ATE for the specified subpopulation.
#'
#' @export
#'
#' @examples
#' # Example data
#' set.seed(123)
#' n <- 100
#' endorse <- sample(0:1, n, replace = TRUE) # Rule compliance (0/1)
#' outcome <- rnorm(n, 50, 10) # Continuous outcomes
#' treat <- sample(0:1, n, replace = TRUE) # Treatment assignment
#' propensity <- runif(n, 0.3, 0.7) # Propensity scores
#'
#' # Calculate ATE for rule compliers
#' rule.cate.neyman(endorse, outcome, treat, group = 'obey')
#'
#' # Calculate ATE for rule non-compliers with propensity weighting
#' rule.cate.neyman(endorse, outcome, treat, group = 'disobey', propensity = propensity)
#'
#' # Calculate ATE for rule compliers with propensity weighting
rule.cate.neyman <- function(endorse.col,outcome,treat,group = 'obey', propensity = NULL){
  if(group =='obey'){ #endorse.col!=0
    index <- which(endorse.col!=0)
  }
  if(group =='disobey'){ #endorse.col==0
    index <- which(endorse.col==0)
  }  
  if(is.null(propensity)) propensity <- rep(1,length(outcome))
  compute.cate.neyman(outcome = outcome[index],
                      treat = treat[index], 
                      propensity = propensity[index])
}

#' Calculate Neyman Estimate of Conditional Average Treatment Effect (CATE)
#'
#' This function computes the Neyman estimate of the Conditional Average Treatment Effect
#' for a subgroup. It can incorporate propensity score weighting to adjust for
#' treatment assignment probabilities when provided.
#'
#' @param outcome A numeric vector of observed outcomes for each unit in the subgroup.
#' @param treat A numeric vector of treatment assignments (typically 1 for treatment, 0 for control)
#' for each unit in the subgroup.
#' @param propensity An optional numeric vector of propensity scores for inverse probability weighting.
#' If NULL (default), no weighting is applied (equivalent to all weights = 1).
#' If provided, should be the same length as outcome and treat.
#'
#' @return A numeric value representing the Neyman estimate of the CATE for the subgroup.
#' This is calculated as the difference between the weighted mean outcome in the
#' treatment group and the weighted mean outcome in the control group.
#'
#' @export
#'
#' @examples
#' # Example without propensity weighting
#' outcome <- c(10, 15, 8, 12, 20, 18, 5, 9, 14, 16)
#' treat <- c(1, 1, 0, 1, 0, 1, 0, 0, 1, 0)
#' compute.cate.neyman(outcome, treat)
#'
#' # Example with propensity weighting
#' propensity_scores <- c(0.6, 0.7, 0.4, 0.8, 0.3, 0.9, 0.2, 0.5, 0.7, 0.4)
compute.cate.neyman <- function(outcome,treat, propensity = NULL){
  if(is.null(propensity)) propensity <- rep(1,length(outcome))
  mean.y1 <- sum(outcome[treat==1]*propensity[treat==1])/sum(treat==1, na.rm = T)
  mean.y0 <- sum(outcome[treat!=1]*propensity[treat!=1])/sum(treat!=1)
  return(mean.y1-mean.y0) 
}

#' Check Covariate Balance Between Treatment and Control Groups
#'
#' This function assesses the balance of covariates between treatment and control groups.
#' It can compute both absolute mean differences (AMD) and standardized mean differences (SMD)
#' for specified variables. Optionally, it can incorporate propensity score weighting.
#'
#' @param data A data.frame containing the dataset with treatment assignment and covariates.
#' @param var.names A character vector specifying the names of the covariates to check for balance.
#' @param propensity An optional numeric vector of propensity scores for weighting.
#' If provided, the covariates will be weighted by these scores before balance checking.
#' Default is NULL (no weighting).
#'
#' @return A data.frame with the following columns for each covariate:
#' \describe{
#' \item{mean.treat}{Mean value in the treatment group (optionally weighted).}
#' \item{mean.control}{Mean value in the control group (optionally weighted).}
#' \item{AMD}{Absolute Mean Difference between treatment and control groups.}
#' \item{10%SMD}{Standardized Mean Difference multiplied by 100 (expressed as a percentage).}
#' }
#'
#' @export
#'
#' @examples
#' # Create example data
#' set.seed(123)
#' n <- 100
#' data <- data.frame(
#' treat = sample(0:1, n, replace = TRUE),
#' age = rnorm(n, 50, 10),
#' income = rnorm(n, 50000, 10000),
#' gender = sample(0:1, n, replace = TRUE)
#' )
#'
#' # Check balance without weighting
#' balance_result <- check.balance(data, c('age', 'income', 'gender'))
#' print(balance_result)
#'
#' # Check balance with propensity score weighting (example)
#' propensity_scores <- runif(n, 0.3, 0.7)
#' balance_weighted <- check.balance(data, c('age', 'income', 'gender'),
#' propensity = propensity_scores)
#' print(balance_weighted)

check.balance <- function(data, var.names, propensity = NULL){
  # data <- data
  if(!is.null(propensity)){data[var.names] <- data[var.names]*propensity}
  X.treat <- data[data$treat == 1,var.names]
  X.control <- data[data$treat != 1,var.names]
  balance <- data.frame(mean.treat =  round(colMeans(X.treat),4),
                        mean.control =  round(colMeans(X.control),4))
  balance$AMD <-  abs(balance$mean.treat-balance$mean.control) #A代表绝对值
  balance['10%SMD'] <- lapply(row.names(balance), function(name){
    var.t <- var(data[data$treat==1,name])
    var.c <- var(data[data$treat!=1,name])
    abs(100*balance[name,'AMD']/sqrt((var.t+var.c)/2))
  }) %>% unlist() 
  return(balance)
}

#' Calculate Empirical Value Function for Treatment Mechanism Evaluation
#'
#' This function computes the empirical value function, which estimates the expected potential outcomes
#' under a specific treatment assignment mechanism. It uses inverse propensity weighting to account
#' for the treatment assignment probabilities.
#'
#' @param reward A numeric vector of observed rewards/outcomes for each unit.
#' @param treat A numeric vector of actual treatment assignments for each unit.
#' @param predictions A numeric vector of predicted treatment assignments (the mechanism to evaluate).
#' @param propensity A numeric vector of propensity scores, i.e., the probability of receiving
#' the actual treatment for each unit. Must be between 0 and 1.
#'
#' @return A numeric value representing the estimated expected potential outcome under the
#' specified treatment assignment mechanism.
#'
#' @export
#'
#' @examples
#' # Example data
#' reward <- c(10, 15, 8, 12, 20) # Observed outcomes
#' treat <- c(1, 0, 1, 0, 1) # Actual treatment assignments
#' predictions <- c(1, 1, 0, 0, 1) # Predicted treatment assignments
#' propensity <- c(0.6, 0.4, 0.7, 0.3, 0.8) # Propensity scores
#'
empirical.value <- function(reward,treat,predictions,propensity){
  numerator <- mean(reward*(predictions==treat)/propensity, na.rm = T)
  denominator <- mean((predictions==treat)/propensity, na.rm = T)
  return(numerator/denominator)
}

#' Compute Error Metrics for CATE Estimates
#'
#' This function calculates various error metrics to evaluate the performance of
#' Conditional Average Treatment Effect (CATE) estimates against the true treatment effects.
#'
#' @param true.tau A numeric vector of the true treatment effects.
#' @param estimate.tau A numeric vector of the estimated treatment effects. Must be the same length as true.tau.
#' @param type A character string specifying the type of error metric to compute.
#' Options are: 'MSE' (Mean Squared Error), 'RMSE' (Root Mean Squared Error),
#' or 'MAE' (Mean Absolute Error). Default is 'MSE'.
#'
#' @return A numeric value representing the computed error metric, rounded to 3 decimal places.
#'
#' @export
#'
#' @examples
#' # Example data
#' true_effects <- c(2.0, 1.5, 3.2, 0.8, 2.5)
#' estimated_effects <- c(1.8, 1.7, 3.0, 1.0, 2.3)
#'
#' # Compute different error metrics
#' compute.Mse.CateEstimates(true_effects, estimated_effects, type = 'MSE')
#' compute.Mse.CateEstimates(true_effects, estimated_effects, type = 'RMSE')
#' compute.Mse.CateEstimates(true_effects, estimated_effects, type = 'MAE')
compute.Mse.CateEstimates <- function(true.tau,estimate.tau, type = 'MSE'){
  if(type == 'MSE') result <- mean((estimate.tau-true.tau)**2)
  if(type == 'RMSE') result <- mean((estimate.tau-true.tau)**2)**0.5
  if(type == 'MAE') result <- abs(estimate.tau-true.tau) %>% mean()
  
  return(round(result, digits = 3))
}


#' Identify Dominant Subgroups for Treatment Efficient Frontier
#'
#' This function identifies dominant subgroups from a set of treatment groups based on
#' their average treatment effect and subgroup size. A subgroup is considered dominant
#' if there exists no other subgroup that has both a strictly larger average treatment effect
#' and a subgroup size at least as large.
#' This is useful for constructing a treatment efficient frontier where we want to
#' identify the most efficient subgroups that cannot be dominated by others.
#'
#' @param average_treatment_effect A numeric vector of average treatment effects for each subgroup.
#' @param subgroup_size A numeric vector of subgroup sizes (number of observations in each subgroup).
#' Must be the same length as average_treatment_effect.
#'
#' @return A numeric vector containing the indices of the dominant subgroups.
#'
#' @export
#'
#' @examples
#' # Example with 5 subgroups
#' ate <- c(2.5, 3.0, 1.8, 3.2, 2.9) # Average treatment effects
#' sizes <- c(100, 80, 120, 70, 90) # Subgroup sizes
#'
#' dominant_indices <- find.dominant.grps(ate, sizes)
#' print(dominant_indices) # Returns indices of dominant subgroups
#'
#' # Display the dominant subgroups' characteristics
#' data.frame(
#' Subgroup = dominant_indices,
#' ATE = ate[dominant_indices],
#' Size = sizes[dominant_indices]
find.dominant.grps <- function(average_treatment_effect,subgroup_size){
  dominant_groups <- numeric()
  for (i in 1:length(average_treatment_effect)){
    is_dominant <- TRUE
    # print(i)
    for (j in 1:length(average_treatment_effect)){
      if (average_treatment_effect[j] > average_treatment_effect[i] & subgroup_size[j] >= subgroup_size[i]){
        is_dominant <- FALSE
        break
      }
    }
    if (is_dominant) {
      dominant_groups <- c(dominant_groups, i)
    }
  }
  return(dominant_groups)
}
