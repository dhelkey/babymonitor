
#' Fit Baby-MONITOR for CPQCC/VON
#'
#' \code{fitBabyMonitor} applies the Baby-MONITOR score  to a single performance indicator. Designed for CPQCC/VON usage.
#'
#' @param minimal_data Data_frame (format must be as specified):
#'
#' 1st column: Outcome vector (0-1 encoding)
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
#' @param num_cat integer. Number of categorical risk adjusters.
#' @param num_cont integer. Number of continuous risk adjusters.
#' @param outcome_type string. Data type of data 'dichotomous'/'cont'. ('cont' requires use_JAGS=TRUE)
#' @param alpha scalar in (0,1). Statistical significance threshold. Posterior (1-alpha)\% intervals are generated.
#' @param use_JAGS bool.  If TRUE, performs inference with statgistical inferance with JAGS (required for continious outcomes)
#' @param prior_var_beta. Prior variance for regression coefficients.
#' @param iters integer. Desired number of posterior MCMC iterations.
#' @param burn_in integer. Number of 'burn-in' MCMC iterations to discard.
#' @param prior_shape.  Prior gamma shape parameter (only for outcome_type=='cont').
#' @param prior_shape.  Prior gamma rate parameter (only for outcome_type=='cont').
#' @param bonferroni logical; if TRUE, posterior intervals are widened with the Bonferroni correction.
#' @param t_scores logical; if TRUE, posterior intervals confidence intervals constructed with the student t-distribution. If FALSE, z-scores are used. Default = TRUE
#' @param dat_out logical; if TRUE, export MCMC iterations and other parameters in the \code{dat} component.
#' @param outcome_na Method for handling any NA values in the outcome vec. 'remove' removes rows with NA outcomes while 'set0' keeps the row and sets the outcome to 0.
#' @param subset_na Method for handling any NA subset values. 'remove' removes rows with NA subset values while 'category' makes a new subset category (coded as 99) for NA values.
#' @param cat_na Method for handling any NA values in the categorical risk adjusters. 'remove' removes rows with NA values while 'category' makes a new category (coded as 99) for NA catorical risk adjusters.
#' @param cont_na Method for handling any NA values in the continous risk adjusters. 'remove' removes rows with NA values while 'median' replaces NA with the median value of the risk adjuster.
#' @param compute_dg logical; if TRUE, compute DG scores, using quantized values of risk adjusters, if applicable.
#'
#' @return
#'	Returns a large list with the following components:
#'
#'		inst_mat: Matrix (one row per institution) computing summary statistics, DGH institution ranking, and intervals with both effect-size and standardized scores.
#'
#' 	     full_subset_mat_baseline, full_subset_mat_nobaseline: Matrix with rankings and intervals for the various subset categories.
#'
#'
#'	    subset_baseline_mat, subset_nobaseline_mat: Matrix with rankings and intervals for subset categories within institution.
#'
#'
#'		group_labels: A vector of the names of each institution
#'
#'	    mcmc_fit: Matrix of MCMC iterations for each coefficient
#'		dg_z: Matrix of computed z score for each institution at each MCMC iterations.
#'
#'		coefs: Names of each coefficient (1st is intercept, then institution, then everything else)
#'
#'		prior_var_vec: Vector of prior variances for each coefficient
#'
#'		model_matrix: Design matrix

fitBabyMonitor = function(minimal_data, num_cat, num_cont,
                          outcome_type = 'dichotomous',
                          alpha = 0.01,
                          use_JAGS = FALSE,
                          subset = FALSE,
                          prior_var_beta = 10,
                          subset_base_catgory = 1,
                          iters = 500,
                          burn_in = 25,
                          dg_gamma = 0.5,
                          prior_shape = 1,
                          prior_rate = 1,
                          bonferroni = TRUE,
                          t_scores = TRUE,
                          outcome_na = 'remove',
                          subset_na = 'category',
                          cat_na = 'category',
                          cont_na = 'median',
                          score_type = 'stat_z',
                          dat_out = FALSE,
                          compute_dg = FALSE)
  {
  if (outcome_type == 'cont'){use_JAGS=TRUE} #Continuous model requires JAGS

  #Process data into standardized format, removing all NA
  dat = parseMinimalData(minimal_data, num_cat, num_cont,
                         subset = subset, outcome_na = outcome_na,
                         subset_na = subset_na, cat_na = cat_na,
                         cont_na = cont_na)

  #Partition data by institution, subset, and institution-subset
  p_inst = partitionSummary(dat$y, dat$inst_vec)
  p_subset = partitionSummary(dat$y, dat$subset)
  p_inst_subset = NULL
  if (subset){p_inst_subset = partitionSummary(dat$y, dat$inst, dat$subset_vec)} #To save time, don't always compute

  #Additive model for risk adjusters: (intercept + categorical variables w/ interactions + continious variables)



  model_mat_cat = modelMatrix(dat$cat_var_mat, interactions = TRUE)
  model_mat_cont = NULL
  if (num_cont > 0){model_mat_cont = modelMatrix(dat$cont_var_mat)}


  model_mat = cbind( rep(1, dat$N), model_mat_cat, model_mat_cont)
  if (num_cat == 0){ 
	model_mat = cbind( rep(1, dat$N), model_mat_cont)
}
  p_tot = dim(model_mat)[2] #Number of parameters in the regression model

  #Store model matrices of different variable types
  dat$model_mat_cat = model_mat_cat
  dat$model_mat_cont = model_mat_cont

  dg_out_vec = NULL
  if (compute_dg){ #Compute D-G score, with categorical variables quantized
    dg_fun = switch(outcome_type,
                    'dichotomous' = designBased,
                    'cont' = designBasedCont)
    dg = dg_fun(dg_gamma, dat$y, dat$pcf_vec_cont, dat$inst_vec)
    dg_out_vec = dg$Z
  }

  #Prior variance vector for regression coefficients
  prior_var_vec = rep(prior_var_beta, p_tot)

  #Fit Model
  ###############################
  if (use_JAGS){#Collect data for JAGS modeling (note that continuous model requires specifying a value for prior_gamma)
    data_list = list(Y = dat$y, X = as.matrix(model_mat), n = dat$N, p = p_tot,
                     prior_var_beta = prior_var_beta)
  }
  if (outcome_type=='dichotomous'){
    if (use_JAGS){
      model_string = "model{
      #Likelihood
      for (i in 1:n){
      Y[i] ~ dbin(prob[i], 1)
      prob[i] = phi(inprod(beta, X[i, ]))
      }
      #Prior
      for (j in 1:p){
      beta[j] ~ dnorm(0, 1/prior_var_beta)
      }
    }"
      #Fit model w/ JAGS
      model_fit = jagsFun(data_list, model_str = model_string, n_iters = iters, n_burnin = burn_in)

      #Extract MCMC iterations from JAGS
      mcmc_iters = as.matrix(model_fit$samples[[1]][ ,1:p_tot])

    } else { #Fit with native R code (quicker)
      mcmc_iters = probitFit(dat$y, model_mat, prior_var_vec,
                             iters = iters + burn_in)[-(1:burn_in),  ]
    }
    #Create Institution level estimates
    #MCMC matrix of individual level probabilities
    p_i_mat = pnorm(as.matrix(tcrossprod(model_mat, mcmc_iters)))

    #Extract posterior row mean, implied observational variance, and row variance
    p_i_vec = apply(p_i_mat, 1, mean)
    p_i_var_vec = apply(p_i_mat, 1, var)
    pq_i_vec = p_i_vec * (1-p_i_vec)
    p_i_overall_var_vec =  pq_i_vec + p_i_var_vec #Law of total variance
  } else if (outcome_type == 'cont'){

    model_string = "model{
      #Likelihood
      for (i in 1:n){
      Y[i] ~ dnorm(mu[i], inv.var)
      mu[i] = inprod(beta, X[i, ])
      }
      #Prior
      inv.var ~ dgamma(prior_shape, prior_rate)
      sigma2 = 1/inv.var
      for (j in 1:p){
      beta[j] ~ dnorm(0, 1/prior_var_beta)
      }
    }"
    data_list$prior_shape = prior_shape
    data_list$prior_rate = prior_rate
    parameters_to_save = c('beta', 'sigma2')
    model_fit = jagsFun(data_list, model_str = model_string, parameters_to_save = parameters_to_save,
                        n_iters = iters, n_burnin = burn_in)

    #Extract MCMC_iters
    mcmc_iters = model_fit$samples[[1]][ ,1:p_tot]
    sigma2_iters = model_fit$samples[[1]][ ,p_tot+1]
    p_i_mat =as.matrix(tcrossprod(as.matrix(model_mat), mcmc_iters))
    p_i_vec = rowMeans(p_i_mat)
    p_i_var_vec = apply(p_i_mat, 1, var)
    sigma2_vec = rep(mean(sigma2_iters), dat$N)
    p_i_overall_var_vec =  sigma2_vec + p_i_var_vec #Law of total variance
  }
  rm(p_i_mat) #Remove to save space

  #Draper-Gittoes-Helkey institution effects
  dgh_inst = dghRank(dat$y, p_i_vec, p_i_overall_var_vec, p_inst$ind_mat)
  dgh_subset = dghRank(dat$y, p_i_vec, p_i_overall_var_vec, p_subset$ind_mat)
  dgh_inst_subset = dghRank(dat$y, p_i_vec, p_i_overall_var_vec, p_inst_subset$ind_mat)

  #Compute score data
  scores_inst = toScore(dgh_inst$D, dgh_inst$S)
  scores_subset = toScore(dgh_subset$D, dgh_subset$S)
  scores_inst_subset = toScore(dgh_inst_subset$D, dgh_inst_subset$S)

  summaryMat = function(part_dat,dgh_mat, score_mat){
    #Combine data into human readable data.frames
    #Choose a score type

    #If null, terminate early
    if (is.null(dgh_mat)){return(NULL)}

    #Extract labels, counts, and observed rate
    out_mat = cbind(part_dat$part_mat,
                    data.frame(n = part_dat$n, o_mean = part_dat$o_mean))

    #Extract effect and SE
    effect_est = dgh_mat$D; effect_se =  dgh_mat$S

    #COmpute Logistic Regression quality ratio estimate
    lr_est = dgh_mat$O / dgh_mat$E

    #Select score est and se for desired standardization method (score_type)
    score_list = switch(score_type,
                        'effect_scaled'= list(score_mat$est_effect_scaled,score_mat$s_effect_scale),
                        'stat_z'= list(score_mat$est_stat_z,score_mat$s_stat_z),
                        'stat_z_scaled'= list(score_mat$est_stat_z_scaled, score_mat$s_stat_z_scaled)
    )
    score_est = score_list[[1]]
    score_se = score_list[[2]]

    out =  cbind(out_mat,
                 data.frame(
                   effect_est = effect_est,
                   effect_se = effect_se,
                   lr_est = lr_est,
                   score_est = score_est,
                   score_se = score_se,
                   stat_z = dgh_mat$Z))
    out[c('dg_z')] = dg_out_vec
    return(out)
  }

  #Summary data
  inst_mat = summaryMat(p_inst, dgh_inst, scores_inst)
  names(inst_mat)[1] = 'inst'

  subset_mat = summaryMat(p_subset, dgh_subset, scores_subset)
  inst_subset_mat = summaryMat(p_inst_subset, dgh_inst_subset,
                               scores_inst_subset)
  if (subset){
    names(subset_mat)[1] = 'subset'
    names(inst_subset_mat)[1:2] = c('inst','subset')
  }

  #Create baseline values w/ toBaseline
  subset_mat_baseline = toBaseline(subset_mat)
  inst_subset_list_baseline = lapply(unique(inst_subset_mat$inst),
                                     function(inst) toBaseline(inst_subset_mat[inst_subset_mat$inst == inst,  ], inst=TRUE))
  inst_subset_mat_baseline = do.call('rbind', inst_subset_list_baseline)

  #Add confidence intervals
  aI = function(mat){#Wrapper for add intervals to include all options
    addIntervals(mat, bonferroni = bonferroni, t_scores = t_scores, alpha = alpha)
  }
  inst_mat = aI(inst_mat)
  subset_mat_nobaseline = aI(subset_mat)
  subset_mat_baseline = aI(subset_mat_baseline)
  inst_subset_mat_nobaseline = aI(inst_subset_mat)
  inst_subset_mat_baseline = aI(inst_subset_mat_baseline)

  if (!dat_out){#Remove saved variables to save storage space.
    dat = list()
  } else{#If dat_out, export variables
    dat$model_mat = model_mat; dat$mcmc_iters = mcmc_iters; dat$prior_var_vec = prior_var_vec
  }
  dat$subset = subset #Save subset info

  return(list(
    dat = dat,
    inst_mat = inst_mat,
    subset_mat_nobaseline = subset_mat_nobaseline,
    subset_mat_baseline = subset_mat_baseline,
    inst_subset_mat_nobaseline = inst_subset_mat_nobaseline,
    inst_subset_mat_baseline = inst_subset_mat_baseline
  ))
}
