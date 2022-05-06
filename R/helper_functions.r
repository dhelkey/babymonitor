#########################
#useful wrapper functions for fitBabyMonitor and jagsUI::jags
babyMonitor = function(minimal_data, num_cat, num_cont,...){
  #' Wrapper function for quick analysis of BabyMONITOR
  #'
  #' \code{linCholSolver} solves Ax = y for x.
  fit = fitBM(minimal_data, num_cat, num_cont,...)
  return(fit$inst_mat)
}

runDG = function(outcome = 'survival', data_use = bm_full,
                 id_var = 'deidhospid',...){
  dgFun = designBased
  outcome_type = 'dichotomous'
  if (outcome=='gv'){
   outcome_type = 'cont'
  }
  #Extract variables from outcome_list
  dat = outcome_list[[outcome]]
  cat_vars = dat$cat_vars
  cont_vars = dat$cont_vars
  num_cat = length(cat_vars)
  num_cont = length(cont_vars)

  #Extract optimal parameter values
  params = outcome_prior_list[[outcome]]
  gamma = params$gamma

  #Extract minimal data
  minimal_data = data_use[ ,c(outcome, id_var, cat_vars, cont_vars)]

  #Draper-Gittoes rankings
  dg = fitDG(minimal_data, num_cat, num_cont,
             gamma = gamma, outcome_type = outcome_type,...)
  return(dg)
}

runBM= function(outcome='survival', data_use = bm_full,
			 id_var = 'deidhospid', ...){
  #' Run BabyMonitor for a given outcome
  outcome_type = 'dichotomous'
  if (outcome=='gv'){
    outcome_type = 'cont'
  }

  #Extract variables from outcome_list
  dat = outcome_list[[outcome]]
  cat_vars = dat$cat_vars
  cont_vars = dat$cont_vars
  num_cat = length(cat_vars)
  num_cont = length(cont_vars)
  minimal_data = data_use[ ,c(outcome, id_var, cat_vars, cont_vars)]
  
  #Extract parameters
  params = outcome_prior_list[[outcome]]
  bm_var_inst = params$inst;
  bm_var_params = params$params

  #Fit & run baby monitor
  bm = fitBM(minimal_data,
					num_cat,
					num_cont,
					outcome_type = outcome_type,
           prior_var_beta_inst = bm_var_inst,
           prior_var_beta_params = bm_var_params,
		   verbose = TRUE, ...)
  return(bm)
}

##############################
computeGreatestChange = function(x1, x2, label = TRUE){
  #' Change in rankings from x1 to x2
  #' if label=FALSE, return raw percent difference

  stopifnot(length(x1)==length(x2))
  pct1 = ecdf(x1)(x1)
  pct2 = ecdf(x2)(x2)
  pct_diff = pct1 - pct2
  #Percentile of the percentile difference
  pct_diff = ecdf(pct_diff)(pct_diff)

  if (!label){
    return(pct_diff)
  }
  cut_vec =  cut(pct_diff, c(0, .1, .9, 1),
                   labels = c('10th Percentile (Drop)',
                              'Middle 80% (~Same)',
                             '90th Percentile (Gain)'))
  return(cut_vec)
}

jagsFun = function(data_list,model_str = '', parameters_to_save = c('beta'),
	n_chains = 1, n_iters = 100, n_burnin = 25, verbose = FALSE){
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
		verbose = verbose)
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

addIntervals = function(data_frame, alpha = 0.01, bonferroni = TRUE,
                        t_scores = TRUE, 
						effect_name = 'effect',
						se_name = 'SE',
						df_name=NULL,
						min_df = 9
						){
  #' Add upper and lower bounds to estimation bounds using
  #' a t-distribtion. If not provided, DF estimated from sample size
  #'
  #'@inheritParams fitBabyMonitor
  if (is.null(data_frame)){return(NULL)}


  if (is.null(df_name)){
	df_vec = data_frame$n - 1
  } else {
	df_vec = data_frame[df_name]
  }

  #Minimum degrees of freedom
  df_vec[df_vec < min_df] = min_df

  #If score_se = 0, the value is fixed, and is not counted in multiple comparisons
  q = sum(data_frame[[se_name]] != 0)
  effect_lower = paste0(effect_name,'_lower')
  effect_upper = paste0(effect_name, '_upper')

  c_values = criticalValues(df_vec, q, alpha = alpha, bonferroni = bonferroni,
                            t_scores = t_scores)


    data_frame[[effect_lower]] = data_frame[[effect_name]] -  
						c_values * data_frame[[se_name]]
						
    data_frame[[effect_upper]] = data_frame[[effect_name]] +
						c_values * data_frame[[se_name]]
  
  data_frame$critical_value = c_values
  return(data_frame)
}

correctInf = function(mat){
  #Convert Inf values to 0
  mat[is.infinite(mat)] = 0
  return(mat)
}

toQuantiles = function(x, derived_levels = 5){
  #' Convert numeric variable to categorical
  #' with derived_levels equally spaced percentiles
  if (length(x) == 0){return(x)}
  q_vec = quantile(as.numeric(x), probs = seq(0, 1, length = derived_levels + 1), na.rm = TRUE)
  q_vec[1] = -Inf
  if (length(unique(q_vec)) < length(q_vec)){
	q_vec = c( sort(unique(floor(q_vec))), Inf)
  }
  labels = 1:(length(q_vec)-1)
  return(as.numeric(as.factor(cut(x, q_vec, labels = labels))))
}
idStr = function(x){ #Convert id (potentially two columns) to vector
		trimws(as.character(paste(x, collapse = '-')))
}
