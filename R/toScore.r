  toScore = function(effect_vec, se_vec){
  #This is just a linear transform, mean and variance transform as well
  #Requires the effect be computed with an estimate and an SE
  if(is.null(effect_vec)){return(NULL)}
	
	#Scaled effects
	sd_effect = sd(effect_vec)
	est_effect_scaled = effect_vec / sd_effect
	s_effect_scaled = se_vec / sd_effect
	
	#Statistical z-scores
	est_stat_z = effect_vec / se_vec
	s_stat_z = se_vec / se_vec #1
		
	#Scaled z-scores
	sd_z = sd(est_stat_z)
	est_stat_z_scaled = effect_vec / (sd_z * se_vec)
	s_stat_z_scaled = se_vec / (sd_z * se_vec)
	
	return(
		data.frame( est_effect_scaled = est_effect_scaled, s_effect_scaled = s_effect_scaled,
	est_stat_z = est_stat_z, s_stat_z = s_stat_z, 
	est_stat_z_scaled = est_stat_z_scaled, 
	s_stat_z_scaled = s_stat_z_scaled)
	)		
  }