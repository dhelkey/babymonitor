#' Fit Baby-MONITOR for CPQCC/VON
#'
#' \code{fitBM} applies the Baby-MONITOR score  to a single performance indicator. Designed for CPQCC/VON usage.
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
#' @param cat_na Method for handling any NA values in the categorical risk adjusters. 'remove' removes rows with NA values while 'category' makes a new category (coded as 99) for NA catorical risk adjusters.
#' @param cont_na Method for handling any NA values in the continous risk adjusters. 'remove' removes rows with NA values while 'median' replaces NA with the median value of the risk adjuster.
#' @return
#'	Returns a large list with the following components:
#'
#'		inst_mat: Matrix (one row per institution)
#'		computing summary statistics, DGH institution ranking,
#'		 and intervals with both effect-size and standardized scores.
#'
#'		group_labels: A vector of the names of each institution
#'
#'	    mcmc_fit: Matrix of MCMC iterations for each coefficient
#'
#'		dg_z: Matrix of computed z score for each institution at each MCMC iterations.
#'
#'		coefs: Names of each coefficient (1st is intercept, then institution, then everything else)
#'
#'		prior_var_vec: Vector of prior variances for each coefficient
#'
#'		model_matrix: Design matrix
fitBM = function(minimal_data, num_cat, num_cont,
                          outcome_type = 'dichotomous',
                          alpha = 0.01,
                          prior_var_beta_inst = 1,
                          prior_var_beta_params = 100,
                          iters = 250,
                          burn_in = 25,
                          prior_shape = 1,
                          prior_rate = 1,
                          bonferroni = TRUE,
                          t_scores = TRUE,
                          outcome_na = 'remove',
                          cat_na = 'category',
                          cont_na = 'median',
                          derived_levels = 5,
                          verbose = FALSE,
						  all_data_out = FALSE)
{
  #Process data into standardized format
  dat = parseMinimalData(minimal_data, num_cat, num_cont,
                        outcome_na = outcome_na,
                         cat_na = cat_na,
                         cont_na = cont_na,
                         derived_levels = derived_levels)

  #Additive statistical model:
  #   (intercept + categorical variables w/ interactions + continious variables)
  intercept_vec = rep(1, dat$N)

  model_mat_cat = modelMatrix(dat$cat_var_mat, interactions = TRUE)

  model_mat_cont = NULL
  if (num_cont > 0){model_mat_cont = modelMatrix(dat$cont_var_mat)}

  #Model matrix containing all risk adjusters + intercept
  model_mat_cov = as.matrix(cbind( intercept_vec, model_mat_cat, model_mat_cont))

  #Number of parameters in the regression model, before institution parameters
  p_tot = dim(model_mat_cov)[2]

  #Model matrix components
  model_mat_inst = as.matrix(t(table(data.frame(hospid = dat$inst_vec,
                                infant = 1:length(dat$inst_vec)))))
  n_inst = dim(model_mat_inst)[2]

  #Assemble matrix defining statistical model
  model_mat =  cbind(model_mat_inst, as.matrix(model_mat_cov))

  prior_var_vec = c(rep(prior_var_beta_inst, n_inst),
                    rep(prior_var_beta_params, p_tot))

  #Set up model and data to fit with JAGS
  model_string = switch(outcome_type,
            'dichotomous'="model{
                        #Likelihood
                        for (i in 1:n){
                          Y[i] ~ dbin(prob[i], 1)
                          prob[i] = phi(inprod(beta, X[i, ]))
                        }
                        #Prior
                        for (j in 1:p){
                          beta[j] ~ dnorm(0, 1/prior_var[j])
                        }
                        }",
        'cont'="model{
            #Likelihood
            for (i in 1:n){
              Y[i] ~ dnorm(mu[i], inv.var)
              mu[i] = inprod(beta, X[i, ])
            }
            #Prior
            inv.var ~ dgamma(prior_shape, prior_rate)
            sigma2 = 1/inv.var
            for (j in 1:p){
              beta[j] ~ dnorm(0, 1/prior_var[j])
            }
}")

  #Specify data in list for JAGS computation
  y = dat$y
  p = dim(model_mat)[2]

  data_list = list(Y = y,
                   X = model_mat,
                   n = length(y),
                   p = p,
                   prior_var = prior_var_vec)
  if (outcome_type=='cont'){
    data_list$prior_shape = prior_shape
    data_list$prior_rate = prior_rate
  }
  
	model_fit = jagsFun(data_list,
					  model_str = model_string,
					  n_iters = iters,
					  n_burnin = burn_in,
					  verbose = verbose)

	#Extract MCMC iterations from JAGS
	mcmc_iters = as.matrix(model_fit$samples[[1]])
	
	#Compute observed statistics from parsed data
	n_vec = as.vector(table(dat$inst_vec))
	O = as.vector(by(y, dat$inst_vec, mean))
	

	#Posterior Samples of Alpha
    inst_iters = mcmc_iters[ ,1:n_inst]
	model_mat_inst = model_mat[ , 1:n_inst]
	
	#Estimate institutional risk from covariates
	cov_iters = mcmc_iters[ ,-(1:n_inst)]
	model_mat_risk =  model_mat[ ,-(1:n_inst)]
	
	#Individual estimated risk
	#Transform function: coefficients -> probability scale
	coefTrans = switch(outcome_type,
		'dichotomous'=function(x) pnorm(x),
		'cont'=function(x) x)

	# MCMC average of
	# (Model matrix w/o institutional effects)* MCMC iterations 
	est_risk_individ = coefTrans(apply(model_mat_risk %*% t(cov_iters),
									   1, mean))
	est_risk_inst = as.vector(by(est_risk_individ, dat$inst_vec, mean))
	
	#Institutional estimates (  E(beta_j | y) / sqrt( V(beta+j | y) )
	D_est = apply(inst_iters, 2, mean)
	D_se = apply(inst_iters, 2, sd)
    Z = D_est / D_se

    inst_mat = data.frame(inst = dat$unique_inst_vec,
						n = n_vec,
                        effect = D_est,
						SE = D_se,
						O=O,
						E=est_risk_inst,
                        Z = Z)	
						
  if (all_data_out){
	return( 
		list(
			inst_mat = inst_mat,
			mcmc_inst = inst_iters,
			mcmc_risk = cov_iters,
			model_mat_inst = model_mat_inst,
			model_mat_risk = model_mat_risk
		)
	)
  } 
  return(inst_mat)
 }
