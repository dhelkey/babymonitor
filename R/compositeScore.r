compositeScore = function(in_list,  score_var = 'Z', se_var = 'z_se'){

n_scores = length(in_list)

n_vars = paste('n',1:n_scores, sep = '')
score_vars = paste(score_var, 1:n_scores, sep='')
se_vars = paste(se_var, 1:n_scores, sep = '')

i=0
for (df in in_list){
  i = i + 1
  df['z_se'] = 1
  names(df)[-1] = paste(names(df)[-1], i, sep='')
  if (i==1){
    out_df=df
  } else {
    out_df = merge(out_df, df, all=TRUE, on='inst')
  }
}
#Number of scores for each institution
#(e.g. some hospitals may be not have data for all quality measures)
n_components = rowSums(!is.na(out_df[n_vars]), na.rm = TRUE)

#Constructing N(0,1) composite score
score_weight_vec = 1/sqrt(n_components)
score_se_vec = 1/n_components

out_df['score'] = rowSums(out_df[score_vars] * score_weight_vec,
                          na.rm=TRUE)

#Compute W-S estimated DF
# Limitations of the Welch-Satterthwaite
# approximation for measurement uncertainty
# calculations
# https://iopscience.iop.org/article/10.1088/0026-1394/37/1/8/
#Equation [1] in Ballico, 2000
DF_num = 1 #By construction (composite sums scores totalling unit variance)
DF_denom = rowSums(score_se_vec^2 * 1/out_df[n_vars], na.rm=TRUE)
out_df['ws'] = DF_num/DF_denom

out_df['total_n'] = rowSums(out_df[n_vars], na.rm = TRUE)
#Regression DF approximation of W-S
out_df['ws_approx'] = out_df$total_n * 0.98

out_df['score_p'] = 1-pt(abs(out_df$score), out_df$ws_approx)
out_df = out_df[order(out_df$inst),]
return(out_df)
}