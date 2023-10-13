#########################################################################################
#                                                                                       #
#                                 DATA GENERATION                                       #
#                                                                                       #
#########################################################################################


set.seed(1)
n = 500 #sample size

betab <- -0.3 # treatment effect
Zb<-rep(c(0,1), each = n/2) # Treatment variable

#Weibull parameters
shape = 0.6 #must be >0
scale = exp(-betab*Zb/shape) #must be >0

#Uniform parameter
theta<-8.543879 # To have ~20% of censored patients

T_tilde <- rweibull(n, shape = shape, scale = scale) # Weibull distribution for event times

C <- runif(n, 0, theta) # Uniform distribution for censoring times

time <- pmin(T_tilde, C)
event <- as.numeric(time ==T_tilde)

simu <- data.frame(patID = 1:n, # patient ID
                   betab = betab,
                   Zb = Zb,
                   time = time, 
                   event = event)
