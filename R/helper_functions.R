#########################
#useful wrapper functions for fitBabyMonitor and jagsUI::jags
babyMonitor = function(minimal_data, num_cat, num_cont,...){
  #' Wrapper function for quick analysis of BabyMONITOR
  #'
  #' \code{linCholSolver} solves Ax = y for x.
  fit = fitBabyMonitor(minimal_data, num_cat, num_cont,
					   compute_dg=TRUE, #Added this
					   ...)
  return(fit$inst_mat)
}

runBM= function(outcome='survival', 
			simple = TRUE,...){
  #' Run BabyMonitor for a given outcome, results in a returner object
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
  bm = fitBabyMonitor(minimal_data, 
					num_cat,
					num_cont,
					outcome_type = outcome_type, 
					use_JAGS = use_JAGS, ...)
  
  if (simple){
	return(bm$inst_mat)
  }
  return(bm)
}

##############################

jagsFun = function(data_list,model_str = '', parameters_to_save = c('beta'), 
	n_chains = 1, n_iters = 100, n_burnin = 25){
	#' Wrapper for jagsUI::jags
	jagsUI::jags(data_list,
		inits = NULL,
		parameters.to.save = parameters_to_save,
		model.file = textConnection(model_str), 
		n.chains = n_chains,
		n.adapt = n_burnin,
		n.iter = n_iters + n_burnin,
		DIC = FALSE, 
		seed = NULL,
		bugs.format = FALSE,
		verbose = TRUE)
}
	
linCholSolver = function(R, y){
  #' Solve System of Equations w/ Cholesky decomposition
  #'
  #' \code{linCholSolver} solves Ax = y for x.
  #' It is first required to take the Cholesky decomposition,
  #' obtaining A = t(R) %*% R. This must be performed outside of this function for speed.
  #' Results should be checked against a known solver, as this function optimizes for speed.
  #'
  #' @param R - Upper triangular matrix of dimension <n x p>; Cholesky decomposition of A.
  #' @param y - Vector of length n.
  q = forwardsolve(t(R), y)
  return(backsolve(R, q))
}

criticalValues = function(n_vec, q = NULL, alpha = 0.01, bonferroni = TRUE, t_scores = TRUE){
  #' Compute critical values for a (100 - alpha)% CI
  #'
  #'@param n_vec Vector of length (q x 1) of sample sizes of the q institutions or categories.
  #' @inheritParams fitBabyMonitor
  if(is.null(q)){q = length(n_vec)}
  adj = 1
  if (bonferroni){adj = q}
  p_interval = 1 - alpha / (2 * adj)
  QT = function(p,n_vec){
    out = rep(0, length(n_vec))
    out[n_vec > 0] = qt(p, n_vec[n_vec > 0])
    return(out)
  }
  if (t_scores){c_values = QT(p_interval, n_vec)
  } else{c_values = rep( qnorm(p_interval), length(n_vec))}
  return(c_values)
}

addIntervals = function(data_frame, alpha = 0.01, bonferroni = TRUE, t_scores = TRUE, composite = FALSE, min_df = 5){
  #' Add upper and lower bounds to score and est bounds...
  #'
  #'@inheritParams fitBabyMonitor
  if (is.null(data_frame)){return(NULL)}

  df_vec = data_frame$n - 1
  df_vec[df_vec == 0] = 1

  #Adjusted degrees of freedom, never less than min_df
  df_vec[df_vec < min_df] = min_df


  if (composite){df_vec = data_frame$ws_approximate_df}


  #If score_se = 0, the value is fixed, and is not counted in multiple comparisons
  q = sum(data_frame$score_se != 0)
  c_values = criticalValues(df_vec, q, alpha = alpha, bonferroni = bonferroni,
                            t_scores = t_scores)

  #Critical values.
  data_frame$score_lower = data_frame$score_est -  c_values * data_frame$score_se
  data_frame$score_upper = data_frame$score_est +  c_values * data_frame$score_se


  if ('effect_est' %in% names(data_frame)){
    data_frame$effect_lower = data_frame$effect_est -  c_values * data_frame$effect_se
    data_frame$effect_upper = data_frame$effect_est +  c_values * data_frame$effect_se
  }
  data_frame$critical_value = c_values

  return(data_frame)

}


correctInf = function(mat){
  #Convert Inf values to 0
  mat[is.infinite(mat)] = 0
  return(mat)
}

toQuantiles = function(x, levels = 4){
  #' Convert numeric variable to categorical 
  #' w/ 4 levels
  if (length(x) == 0){return(x)}
  q_vec = quantile(as.numeric(x), probs = seq(0, 1, length = levels + 1), na.rm = TRUE)
  q_vec[1] = -Inf
  if (length(unique(q_vec)) < length(q_vec)){
	q_vec = c( sort(unique(floor(q_vec))), Inf)
  }
  labels = 1:(length(q_vec)-1)
  return(as.numeric(as.factor(cut(x, q_vec, labels = labels))))
}

idStr = function(x){ #Make sure id (potentially two columns) is a vector
  as.character(paste(x, collapse = '-'))
}












































