bayesianRegression = function(y, X,  prior_var_vec, iters = 100){
	#' Generate MCMC Samples of Bayesian Linear Regression
	#'
	#' \code{bayesianRegression} generates MCMC samples from Bayesian Linear Regression.
	#' Assumes prior mean centered at 0.
    #' Written by Daniel Kirsner
	#' @param y <nX1> observed data vector
	#' @param X <nXp> design matrix (takes objects of class matrix or sparce Matrix)
	#' @param prior_var_vec <px1> prior variance vector
	#' @param iters Integer; number of desired output iterations
    n = length(y)
	p = dim(X)[2]  	
	stopifnot(length(prior_var_vec) == p)
	
	#Data augmentation
	X_aug = as.matrix(rbind(X, diag(p)))
	y_aug = c(y, rep(0, p))
	y_aug = as.numeric(y_aug)
	
	bayesian_model=lm(y_aug~0+X_aug) #, weights = 1/prior_var_vec) 
	mcmc_mat=MASS::mvrnorm(iters,coef(bayesian_model), vcov(bayesian_model))
    return(mcmc_mat)
}