#########################################################################################
#                                                                                       #
#                              DATA ANALYSIS - COX                                      #
#                                                                                       #
#########################################################################################



### COX model ###
res <- summary(survival::coxph(formula = survival::Surv(time, event) ~Zb, data = simu))

COX_beta_hat <- res$coefficients[1,1];COX_beta_hat
COX_se_hat <- res$coefficients[1,3]; COX_se_hat