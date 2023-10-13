#########################################################################################
#                                                                                       #
#                              DATA ANALYSIS - GEE                                      #
#                                                                                       #
#########################################################################################

### load packages ###
library(pseudo)
library(geepack)


### Defining the time points used to compute pseudo-observations ###
K = 5 # Number of time points


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

b$ipseudo <- 1-b$pseudo # needed because cloglog in geese is defined as log(-log(1-y))
b$event<- as.factor(b$event)
  
### Fit a GEE model ###
fit_GEE <- geepack::geese(ipseudo ~ as.factor(Zb) + as.factor(tpseudo),
                          id = id, data = b, scale.fix=TRUE, family=gaussian, jack = TRUE,
                          mean.link = "cloglog", corstr="independence")

as.data.frame(cbind(mean = fit_GEE$beta, 
                    SE_ajust = sqrt(diag(fit_GEE$vbeta.ajs)), 
                    SE_sandwich = sqrt(diag(fit_GEE$vbeta))))
