#########################################################################################
#                                                                                       #
#                           DATA ANALYSIS - Bayesian GMM                                #
#                                                                                       #
#########################################################################################

### load packages ###
library(tidyverse)
library(pseudo)
library(rstan)



### Defining the time points used to compute pseudo-observations ###
K = 5 # Number of time points
np = K+1
n = length(simu$patID)

event_times <- simu$time[simu$event == 1]
cutoffs <- quantile(event_times, seq(0, 1, length.out = K+2)[2:(K+1)]) # take the quantile of the event time

### Compute the pseudo-observations ###
pseudo <- pseudo::pseudosurv(simu$time, simu$event, tmax=cutoffs)

#arrange the data
b <- NULL
for(it in 1:length(pseudo$time)){
  b <- rbind(b,cbind(simu,pseudo = pseudo$pseudo[,it],
                     tpseudo = pseudo$time[it],id=1:nrow(simu)))
}
b <- b[order(b $id),]
b$event<- as.factor(b$event)
  

  
b$tpseudo1 <- 1*(b$tpseudo==pseudo$time[[1]])
b$tpseudo2 <- 1*(b$tpseudo==pseudo$time[[2]])
b$tpseudo3 <- 1*(b$tpseudo==pseudo$time[[3]])
b$tpseudo4 <- 1*(b$tpseudo==pseudo$time[[4]])
b$tpseudo5 <- 1*(b$tpseudo==pseudo$time[[5]])

  
  
 
X = matrix(c(rep(1,n*K), rep(simu$Zb, each = K), tpseudo2 = b$tpseudo2, tpseudo3 = b$tpseudo3,tspeudo4 = b$tpseudo4, tpseudo5 = b$tpseudo5),
             nrow = n*K, ncol = np)
data <- list(X = X, Y = b$pseudo, n = n, N = n, K = K, np = np)


### fit a Bayesian GMM with independence working matrix ###
GMM_ind <- rstan::stan_model("Bayesian_GMM_on_simulated_pseudo_vector.stan") # compile Stan model
  
# Calculate starting parameters for the 3 chains
b <- b%>%
  mutate(pseudo_cut1 = 0.99*(pseudo>=0.99) + 0.01*(pseudo<=0.01)+ pseudo*(pseudo<0.99 & pseudo>0.01),
         pseudo_cut2 = 0.95*(pseudo>=0.95) + 0.05*(pseudo<=0.05)+ pseudo*(pseudo<0.95 & pseudo>0.05),
         pseudo_cut3 = 0.90*(pseudo>=0.90) + 0.10*(pseudo<=0.10)+ pseudo*(pseudo<0.90 & pseudo>0.10))
  
fit_LM_cloglog1 <- lm(log(-log(pseudo_cut1)) ~ as.factor(Zb)+as.factor(tpseudo), data = b)
fit_LM_cloglog2 <- lm(log(-log(pseudo_cut2)) ~ as.factor(Zb)+as.factor(tpseudo), data = b)
fit_LM_cloglog3 <- lm(log(-log(pseudo_cut3)) ~ as.factor(Zb)+as.factor(tpseudo), data = b)
  
fit_BayesGMM_ind <- rstan::sampling(GMM_ind, data = data, chains = 3, iter = 6000, warmup = 1000, thin = 5, cores = 3,
                                    init = list(chain1 = list(beta = as.numeric(fit_LM_cloglog1$coef)),
                                                chain2 = list(beta = as.numeric(fit_LM_cloglog2$coef)),
                                                chain3 = list(beta = as.numeric(fit_LM_cloglog3$coef))),
                                    save_warmup =T, seed = 1)
summary(fit_BayesGMM_ind)$summary


### fit a Bayesian GMM with exchangeable working matrix ###
GMM_exch <- rstan::stan_model("Bayesian_GMM_on_simulated_pseudo_exch_wcm_vector.stan") # compile Stan model

fit_BayesGMM_exch <- rstan::sampling(GMM_exch, data = data, chains = 3, iter = 6000, warmup = 1000, thin = 5, cores = 3,
                                     init = list(chain1 = list(beta = as.numeric(fit_LM_cloglog1$coef)),
                                                 chain2 = list(beta = as.numeric(fit_LM_cloglog2$coef)),
                                                 chain3 = list(beta = as.numeric(fit_LM_cloglog3$coef))),
                                     save_warmup =T, seed = 1)
summary(fit_BayesGMM_exch)$summary

### fit a Bayesian GMM with exchangeable working matrix ###
GMM_ar1 <- rstan::stan_model("Bayesian_GMM_on_simulated_pseudo_ar1_wcm_vector.stan") # compile Stan model

fit_BayesGMM_ar1 <- rstan::sampling(GMM_ar1, data = data, chains = 3, iter = 6000, warmup = 1000, thin = 5, cores = 3,
                                     init = list(chain1 = list(beta = as.numeric(fit_LM_cloglog1$coef)),
                                                 chain2 = list(beta = as.numeric(fit_LM_cloglog2$coef)),
                                                 chain3 = list(beta = as.numeric(fit_LM_cloglog3$coef))),
                                     save_warmup =T, seed = 1)
summary(fit_BayesGMM_ar1)$summary
