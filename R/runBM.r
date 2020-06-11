runBM= function(outcome='survival', inst_mat = TRUE){
  #' Run BabyMonitor for a given outcome, results in a returner object
  #' Don't upload to Github
  outcome_type = 'dichotomous'
  use_JAGS=FALSE
  if (outcome=='gv'){
    outcome_type = 'cont'
	use_JAGS=TRUE
  }
  dat = outcome_list[[outcome]]
  cat_vars = dat$cat_vars
  cont_vars = dat$cont_vars
  num_cat = length(cat_vars)
  num_cont = length(cont_vars);
  id_var = 'hospid'
  minimal_data = neonatal[ ,c(outcome, id_var, cat_vars, cont_vars)]

  #Fit baby monitor
  bm = fitBabyMonitor(minimal_data, num_cat, num_cont, outcome_type = outcome_type, use_JAGS = use_JAGS)
  
  if (inst_mat){
	return(bm$inst_mat)
  }
  return(bm)
}


