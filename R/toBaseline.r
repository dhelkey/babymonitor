toBaseline = function(data_frame, inst = FALSE){
#' Should be called before generating intervals
	return_cols = c('subset', 'n', 'o_mean', 'effect_est', 'effect_se', 'score_est', 'score_se', 'stat_z')
	if (inst){return_cols = c('inst', return_cols)}

	if (is.null(data_frame)){return(NULL)}
	
	#Identify if the first subset is present
	base_index = which(data_frame$subset == 1)
	if (length(base_index) == 0){ #First subset not present
		df = data.frame(matrix(0, nrow = 0, ncol = length(return_cols)))
		names(df) = return_cols
		return(df)
	}
	#If it is present, then change the relevent columns
	sub_effect = data_frame$effect_est[base_index]
	sub_score = data_frame$score_est[base_index]
	
	data_frame$effect_est = data_frame$effect_est - sub_effect

	data_frame$score_est = data_frame$score_est - sub_score

	data_frame$score_se[base_index] = data_frame$effect_se[base_index] = 
		data_frame$stat_z[base_index] = 0

	return(data_frame[ ,return_cols])
}
  