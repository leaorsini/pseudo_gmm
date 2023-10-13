#########################################################################################
#                                                                                       #
#                 DATA ANALYSIS - Bayesian Piecewise exponential model                  #
#                                                                                       #
#########################################################################################

### load packages ###
library(spBayesSurv)


##Interval Murray :
r = sum(as.numeric(simu$event)) #event number
M = floor(max(5, min(r/8, 20)))


res <- summary(spBayesSurv::indeptCoxph(formula = survival::Surv(time, event) ~Zb, data = simu,
                                        prior=list(M=M, r0=1, beta0 = 0, S0 = 10000),
                                        mcmc=list(nburn=2000, 
                                                  nsave=8000, 
                                                  nskip=0, 
                                                  ndisplay=1000)))
PEM_beta_hat <- res$coef[1,1];PEM_beta_hat
PEM_se_hat <- res$coef[1,3]; PEM_se_hat
