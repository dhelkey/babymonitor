#' Parse input data to fitBM
#'
#' @inheritParams fitBabyMonitor
parseLite = function(minimal_data, num_cat, num_cont,
                            outcome_na = 'remove',
                            cat_na = 'category',
                            cont_na = 'median',
                            unknown_category_code=99){

  minimal_data = as.matrix(minimal_data)

  #Count total records
  N_full = dim(minimal_data)[1]

  #Sort by institution
  minimal_data = minimal_data[order(minimal_data[ ,2]), ]
  if (sum(!complete.cases(minimal_data[ ,2])) > 0){
    stop('Missing values for institution not allowed')}

  #NA outcome variables
  if (outcome_na == 'remove'){ #Remove by default
    minimal_data = minimal_data[complete.cases(minimal_data[ ,1]), ]
  } else if (outcome_na == 'set0'){ #Set to 0
    minimal_data[!complete.cases(minimal_data[ ,1]),1] = 0
  }

  #Categorical risk adjusters
  cat_var_mat = cat_var_locat = NULL
  if (num_cat > 0){
    cat_var_locat = 3:(3 + num_cat - 1)
    #Remove NA categoricals
    if (cat_na == 'remove'){
      minimal_data = minimal_data[
        complete.cases(minimal_data[ ,cat_var_locat]), ]
    }
  }

  #Continuous risk adjusters (set missing values to median)
  imputeFun = function(x){
    if (is.numeric(x)){
      x[!complete.cases(x)] = median(x, na.rm = TRUE)
    }
    return(x)
  }

  cont_var_mat = NULL
  if (num_cont > 0){
    cont_var_locat = (3 + num_cat):(3 + num_cat + num_cont - 1)
    if (cont_na == 'remove'){
      minimal_data = minimal_data[
        complete.cases(minimal_data[ ,cont_var_locat]), ]
    }
    cont_var_mat = as.data.frame(minimal_data[  ,cont_var_locat])
    if (cont_na == 'median'){
      cont_var_mat = cbind(sapply(cont_var_mat, imputeFun))
    }
    colnames(cont_var_mat) = names(minimal_data)[cont_var_locat];
  }
  class(cont_var_mat) <- "numeric"

  #This is done last, after we deal with all NA values
  #Extract categoricals as a matrix and explicitly turn into a factor
  if (num_cat > 0){
    cat_var_mat =as.matrix(minimal_data[  ,cat_var_locat])
    if (cat_na == 'category'){
      cat_var_mat[is.na(cat_var_mat)] = toString(unknown_category_code)
    }
    colnames(cat_var_mat) = names(minimal_data)[cat_var_locat]

    #Explicitly turn each categorical into a factor
    for (i in 1:num_cat){
      cat_var_mat[ ,i] = as.factor(cat_var_mat[ ,i] )
    }
  }

  #Extract variables
  inst_vec = as.factor(minimal_data[ ,2])
  y = as.numeric(minimal_data[ ,1])

  return(list(
    indicator_name = names(minimal_data)[1],
    o_overall = mean(y),
    y = y,
    N_full = N_full,
    N = length(y),
    p = length(unique(inst_vec)),
    inst_vec = inst_vec,
    cat_var_mat = cat_var_mat,
    cont_var_mat = cont_var_mat,
    unique_inst_vec = unique(minimal_data[ ,2])))
}
