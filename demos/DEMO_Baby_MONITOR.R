#Run Baby-Monitor and Draper-Gittoes manually
# fitBM() and fitDG()

#Available datasets:
# bm_full: infants with birth year 2008-2018
# bm_train:     ""                 2013-2015
# bm_test:    ""                   2016-2018

#Uncomment to install the latest babymonitor code from github
# devtools::install_github('dhelkey/babymonitor', force = TRUE)
# library('babymonitor')


outcome = 'survival'
iters = 1500
par(mfrow = c(2,2))
id_var = 'deidhospid'
data_use = bm_test #

outcome_type = 'dichotomous'
if (outcome=='gv'){
  outcome_type = 'cont'
}

#Extract variables from outcome_list
dat = outcome_list[[outcome]]
cat_vars = dat$cat_vars
cont_vars = dat$cont_vars
print(cat_vars)
print(cont_vars)
num_cat = length(cat_vars)
num_cont = length(cont_vars)
minimal_data = data_use[ ,c(outcome, id_var, cat_vars, cont_vars)]

#Tuning parameters (saved in outcome_prior_list)
params = outcome_prior_list[[outcome]]
bm_var_inst = params$inst; #Prior varance - institutional random effect
bm_var_params = params$params #Prior variance -
dg_gamma = params$gamma

#Fit & run baby monitor, including option to return all data (including MCMC iterations)
bm_fit = fitBM(minimal_data,
               num_cat,
               num_cont,
               outcome_type = outcome_type, verbose = TRUE,
               prior_var_beta_inst = bm_var_inst,
               prior_var_beta_params = bm_var_params,
               iters=iters,
               all_data_out=TRUE)
#TEST
par(mfrow = c(2,2))
bm = bm_fit$inst_mat
hist(bm$O)
hist(bm$Z)
hist(bm$effect)
hist(bm$E)

#Compare w/ Draper-Gittoes
dg = fitDG(minimal_data, num_cat, num_cont,
           gamma = dg_gamma, outcome_type = outcome_type)

#Diagnostic plots - 122 NICUs
par(mfrow = c(2,3))
plot(bm$n, bm$effect,
     main = 'NICU size vs estimated effect',
     xlab  = 'Total infants',
     ylab = 'Estimated Effect')


plot(bm$n, bm$Z,
     main = 'Institution size vs Z',
     xlab  = 'Total infants',
     ylab = 'Baby-MONITOR Z')


plot(bm$O, bm$E,
     main = paste('Observed vs Expected - \n', outcome),
     xlab = 'Observed Percentage',
     ylab = 'Expected Percentage')
abline(0, 1, col = 'red', lwd = 2)

plot(bm$O, bm$Z,
     main = 'Observed vs Z',
     xlab = 'Observed',
     ylab = 'Z')
abline(0, 1, col = 'red', lwd = 2)

plot(dg$Z, bm$Z,
     xlab = 'Draper-Gittoes Z',
     ylab = 'Baby-MONITOR Z',
     main = 'Draper Gittoes vs Baby-Monitor Z-Scores')
abline(0, 1, col = 'red', lwd = 2)

#Extracting MCMC iterations (if desired)
#Requires setting all_data_out=TRUE in fitBM()

#MCMC iterations:
#Institutional effects (alpha_i's)
#mcmc_inst: n_institutions x MCMC iterations
#model_mat_inst: Indicator matrix of individual institutional membership
mcmc_inst = bm_fit$mcmc_inst
model_mat_inst = bm_fit$model_mat_inst
print(dim(mcmc_inst))

#MCMC iterations and model matrix for all other parameters (including intersept)
mcmc_risk = bm_fit$mcmc_risk
model_mat_risk = bm_fit$model_mat_risk


##########################################################
#Evaluate outcomes using provided wrapper functions:
# runBM(), runDG()
par(mfrow = c(3,3))
dg_list = list()
bm_list = list()
for (outcome in outcomes){
  print( paste('Computing - ', outcome))
  dg = runDG(outcome, data_use = data_use)
  bm = runBM(outcome, iters = iters, data_use = data_use)
  dg_list[[outcome]] = dg
  bm_list[[outcome]] = bm
  plot(dg$Z, bm$Z, main=outcome )
  abline(0,1, col = 'red')
}


#Construct composite score using each methodology
# bm - Bayesian Baby-MONITOR methodology
# dg - Draper-Gittoes, quantizing continious variables (5 quantiles)
c_bm = compositeScore(bm_list)
c_dg = compositeScore(dg_list)


#Uncomment to load saved composite data
#(computing Baby-MONITOR is time intensive)
#c_bm = c_bm_saved
#c_dg = c_dg_saved


#Viewing important variabes (score, p-value based on t-distributed scores,
# score components 1-9)
head(c_bm[order(c_bm$score_p), c('score', 'score_p', 'Z1','Z2',  'Z3','Z4','Z5','Z6','Z7','Z8','Z9')])

par(mfrow = c(1,1))
plot(c_bm$score, c_dg$score,  xlab = 'Baby-MONITOR Composite',
     ylab = 'Draper-Gittoes composite')
abline(0,1, col = 'red')

#Appendix
###############################################

#Composite score analysis

# Welch–Satterthwaite approxomite degrees of freedom
#Total number of records
#(e.g. if infant has a record for ante and a record for survival - 2 records) vs
#estimated DF from W-S
plot(c_bm$total_n, c_bm$ws) #Note the outlier
abline(0,1, col = 'red')
#Identifying problem (one estimated)
prob_inst = c_bm[c_bm$total_n >2000 & c_bm$ws <300, ]
#1 institution only has 1 record for antenatal steroids
#table(data_use[data_use$deidhospid==prob_inst$inst, c('ante') ])
#Unbalanced sample sizes known to cause problems (Ballico, 2000))
#https://iopscience.iop.org/article/10.1088/0026-1394/37/1/8/
n_vars = paste('n',1:9, sep = '')
n_vars_use = n_vars[n_vars %in% names(prob_inst)]
prob_inst[ ,n_vars_use]


#Regression to identify approximate W-S (98% of "total records")
#  after removing outlier
ws_data = data.frame(y = c_bm$ws, x = c_bm$total_n)
outlier_index = (ws_data$y < 400) & (ws_data$x > 2000)
lm('y ~ x -1', data = ws_data[!outlier_index, ] )


