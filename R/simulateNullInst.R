simulateNullInst = function(minimal_data,
							 num_cat, num_cont, var_intercept = 40,
                             var_cat = 10,
                             var_cat_interaction = 10, var_cont = 10, 
							 iters = 100){
  #Uses the same default parameters as fitBabyMonitor


  ###Run through fitBabyMonitor to parse...thereby using the same default options :)
  r = fitBabyMonitor(minimal_data, num_cat, num_cont,
                     var_intercept = var_intercept, var_cat = var_cat,
                     var_cat_interaction = var_cat_interaction,
                     var_cont = var_cont,
                     dat_out = TRUE,
                     n_cutoff = 1)

  #Extract variables
  N=r$dat$N #total number of individuals
  p = r$dat$p #total number of institutions
  inst_mat = r$inst_mat[ ,c('inst', 'n', 'o_mean')]

  #Posterior mean from fit data
  beta_hat = colMeans(as.matrix(r$dat$mcmc_iters))
  ind_probit_vec = pnorm(as.numeric(r$dat$model_mat %*% beta_hat))

  #Run null simulations in a loop, generating data at each point
  stat_z_mat = effect_mat = effect_lower_mat = effect_upper_mat =
  effect_se_mat = score_mat = matrix(NA, nrow = p, ncol = iters)
  
  data_mat = matrix(NA, nrow = N, ncol = iters)
  
  for (i in 1:iters){#Loop over monte carlo iterations
  
	#Simulate data
    data_mat[ ,i] = rbinom(N, 1, ind_probit_vec)

    m_temp = minimal_data
    m_temp[ ,1] = data_mat[ ,i]
    r_temp =  fitBabyMonitor(m_temp, num_cat, num_cont, var_intercept = var_intercept,
                             var_cat = var_cat,
                             var_cat_interaction = var_cat_interaction, var_cont = var_cont,
                             n_cutoff = 1,
							  bonferroni = FALSE)

    stat_z_mat[ ,i] = r_temp$inst_mat$stat_z
    effect_mat[ ,i] = r_temp$inst_mat$effect_est
	effect_lower_mat[ ,i] = r_temp$inst_mat$effect_lower
	effect_upper_mat[ ,i] = r_temp$inst_mat$effect_upper
    effect_se_mat[ ,i] = r_temp$inst_mat$effect_se
    score_mat[ ,i] = r_temp$inst_mat$score_est
  }

  return( list(
    inst_mat = inst_mat,
    stat_z = stat_z_mat,
    effect = effect_mat,
	effect_lower = effect_lower_mat,
	effect_upper = effect_upper_mat,
    effect_se = effect_se_mat,
    score = score_mat,
	data_mat = data_mat
  ))

}
