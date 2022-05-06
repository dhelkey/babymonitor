#' Draper-Gittoes scores for dataframe input
#'
#' \code{fitDG} applies the Draper-Gittoes score to a single performance indicator. Designed for CPQCC/VON usage.
#'
#' @param minimal_data Data_frame (format must be as specified):
#'
#' 1st column: Outcome vector (0-1 or continious)
#'                 (must set outcome_type = 'cont' for continious outcomes)
#'
#' 2nd column: Institution ID
#'
#' 3rd column: Subset (if subset == TRUE):
#'
#' Next: num_cat columns of categorical variables (num_cat can equal 0)
#'
#' Next: num_cont columns of continuous variables (num_cont can equal 0)
#'
#' @param subset logical; if TRUE perform analysis with a subset variable ( and \code{minimal_data} must have subset data in the 3rd column)
#'
fitDG = function(minimal_data,
                 num_cat, num_cont,
                 subset = FALSE,
                 outcome_type = 'dichotomous',
                 gamma = 0.5,
                 outcome_na = 'remove',
                 subset_na = 'category',
                 cat_na = 'category',
                 cont_na = 'median',
                 derived_levels = 5
                 ){
 #Parse w/ passed options
  dat = parseMinimalData(minimal_data, num_cat, num_cont,
                         subset = subset,
                         outcome_na = outcome_na,
                         subset_na = subset_na,
                         cat_na = cat_na,
                         cont_na = cont_na,
                         derived_levels = derived_levels)

  dgFun = designBased
  if (outcome_type=='cont'){
    dgFun = designBasedCont
  }
  
  dg = dgFun(gamma, dat$y, dat$pcf_vec_cont,dat$inst_vec)
  
  return(dg)
}

###Test
# minimal_data = bm_full[ ,c('survival','deidhospid','male', 'cmal', 'gaweeks')]
# t = runDG(minimal_data, 2, 1)
