modelMatrix = function(x,
  interactions = FALSE, intercept = FALSE){
	#' Create a sparse design matrix from data frame
	#'
	#' \code{modelMatrix} takes in a data.frame and encodes it as a sparse model.Matrix object.
	#' Factors use treatment ("one-hot") encoding, creating indicator variables from categorical
	#' and binary variables.
    #' Note that names of categorical variables will change slightly
	#'	(n-1 unique values appended to name)
	#'
	#' @param x  A data.frame with finite, non-missing elements.
	#' @param interactions Logical; if \code{TRUE} then all interactions included.
	#' @param intercept Logcal; if \code{TRUE} then include a leading column of 1's.
	if(is.null(x))return(NULL) #Dealing w/ empty data
	#if (sum(sapply(x, is.finite)) != 0){
	#	stop('Creating model matrix requires all input elements be finite')}

	#Set R to use 1-hot encoding
	options(contrasts = rep("contr.treatment", 2))
    model_formula = as.formula( ~.) #All first order Effects
    if (interactions){model_formula = as.formula( ~(.)^2)} #1st order and interactions
   mat = model.matrix(model_formula, data = data.frame(x))

    if (!intercept){
		#Remove intercept column, making sure names are preserved
		n_vec = colnames(mat)
		mat = data.frame(mat[ ,-1])
		colnames(mat) = n_vec[-1]
	}
    #Drop any empty levels
	n_vec = colnames(mat)
	nonempty_indices = apply(abs(mat),2, sum) != 0
	mat = data.frame(mat[ ,nonempty_indices])
	p = length(mat)
	colnames(mat) = n_vec[nonempty_indices]
	n_cols = length(colnames(mat))

	#Export as Matrix
	mat = Matrix(sapply(mat, as.numeric))
   return( mat)
}
