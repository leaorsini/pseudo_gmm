#########################################################################################
#                                                                                       #
#                         DATA ANALYSIS - Frequentist GMM                               #
#                                                                                       #
#########################################################################################

### load packages ###
library(pseudo)
library(MASS)

### R function to implement the freq GMM with a cloglog link function ###
gmm_cloglog <- function (formula = formula(data), id = id, data = parent.frame(), b = NULL,
                         tol = 1e-8, maxiter = 1000, family = gaussian, corstr = "independence"){
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$b <- m$tol <- m$maxiter <- m$link <- m$varfun <- m$corstr <- m$family <- m$invfun <- NULL
  
  if (is.null(m$id))
    m$id <- as.name("id")
  
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  y <- as.matrix(model.extract(m, "response"))
  x <- model.matrix(Terms, m)
  
  # get class id
  id <- model.extract(m, id)
  
  # get number of observe persons
  nobs <- nrow(x)
  
  # get number of parameters, include (intercept)
  np <- ncol(x)
  
  # get the name of X matrix colomn (variables)
  xnames <- dimnames(x)[[2]]
  
  if (is.character(family))
    family <- get(family)
  if (is.function(family))
    family <- family()
  
  # Starting parameters: if b is not NULL then sign b to beta; otherwise, use gml to get initial beta
  if (!is.null(b)) {
    beta <- matrix(as.double(b), ncol = 1)
    if (nrow(beta) != np) {
      stop("Dim beta != ncol(x)")
    }
  }
  else {
    # message("\n","running glm to get initial regression estimate")
    mm <- match.call(expand.dots = FALSE)
    mm$b <- mm$tol <- mm$maxiter <- mm$link <- mm$varfun <- mm$corstr <- mm$id <- mm$invfun <- NULL
    mm[[1]] <- as.name("glm")
    beta <- eval(mm, parent.frame())$coef
    
    beta <- as.numeric(beta)
  }
  if (length(id) != length(y))
    stop("Id and y not same length")
  
  #get the maximum number of iteration
  maxiter <- as.integer(maxiter)
  
  # Correlation structure variable CORR must be one of IND, EXCH or AR-1
  corstrs <- c("independence", "exchangeable", "AR-1")
  
  dist <- family$family
  
  ### start sign value to calculation ###
  #### define:
  # x- X matrix,
  # y- response matrix
  # beta - coefficients of X variables
  # id- id list
  # uid- unique id
  # nobs- number of rep time
  # nsub- number of unique id
  # np - number of parameter
  # m0 - m0 matrix, nobs*nobs matrix, depend on correlation structure
  # m1 - m1 matrix, nobs*nobs matrix, depend on correlation structure
  
  # U - np*1 matrix 
  # C - np*np matrix
  # U - np*1 matrix
  # C - np*np matrix
  # ui - np*1 matrix for i-th id
  # Ufirstdev - np*np matrix for first derivative of U
  # firstdev - np*np matrix of first derivative 
  
  # for each id:
  # xi- i-th x
  # yi - i-th y
  # ni- number of row of xi
  
  # mui- i-th mu
  # mui_dev - derivate of mui
  # vui - marginal variance for i-th id
  
  ### Calculation
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  obs <- lapply(split(id, id), "length")
  nobs <- as.numeric(obs)
  nsub <- length(nobs)
  np <- dim(x)[[2]]
  
  ################ 	GMM iteration #######
  betadiff <- 1
  iteration <- 0
  betanew <- beta
  
  # GMM iteration
  while(betadiff > tol && iteration < maxiter)
  {
    
    # initial value
    beta <- betanew
    if (corstr == "independence") {
      U <- matrix(rep(0,np),nrow=np)
      C <- matrix(rep(0,np*np),nrow=np)
      ui <- matrix(rep(0,np),nrow=np)
      Ufirstdev <- matrix(rep(0,np*np),nrow=np)
      firstdev <- matrix(rep(0,np*np),nrow=np)
    }
    else {
      U <- matrix(rep(0,2*np),nrow=2*np)
      C <- matrix(rep(0,2*np*2*np),nrow=2*np)
      ui <- matrix(rep(0,2*np),nrow=2*np)
      Ufirstdev <- matrix(rep(0,2*np*np),nrow=2*np)
      firstdev <- matrix(rep(0,2*np*np),nrow=2*np)
    }
    # one iteration
    
    
    
    
    # loc1 - start location for row in xi
    # loc2 - end location for row in xi
    loc1 <- 0
    loc2 <- 0
    for (i in 1:nsub)
    {
      # set start location for next xi
      loc1 <- loc2+1
      loc2 <- loc1+nobs[i]-1
      
      yi <- as.matrix(y[loc1:loc2,])
      xi <- x[loc1:loc2,]
      ni <- nrow(yi)
      
      # set m0, m1
      m0 <- diag(ni)
      # set m1 by corr structure
      if (corstr == "independence") {
        m1 <- matrix(rep(0,ni*ni),ni)
      }
      else if (corstr == "exchangeable") {
        m1 <- matrix(rep(1,ni*ni),ni) - m0
      }
      else if (corstr == "AR-1") {
        m1 <- matrix(rep(0,ni*ni),ni)
        for (k in 1:ni) {
          for (l in 1:ni) {
            if (abs(k-l)==1) m1[k,l] <-1
          }
        }
      }
      
      # change ui, mui_dev, vui depending on distribution
      if (dist == "gaussian") {
        mui <- exp(-exp(xi %*% beta)) # using a cloglog link function 
        mui_dev <- diag(as.vector(-exp(xi%*%beta)))%*%diag(as.vector(mui))  
        vui <- diag(ni)     # identity matrix for marginal variance (gaussian case)
      }
      
      # calculate mui, wi, zi, c, C, U, di, firstdev, Ufirstdev
      # depending on corr structure
      if (corstr == "independence") {
        wi <- t(xi) %*% mui_dev %*% vui %*% m0 %*% vui
        ui0 <- (1/nsub)*wi %*% (yi-mui)
        ui[1:np,] <- ui0
        C <- C + ui %*% t(ui)
        U <- U + ui
        
        di0 <- -(1/nsub) * wi %*% mui_dev %*% xi
        firstdev[1:np,] <- di0
        Ufirstdev <- Ufirstdev + firstdev
      }
      
      
      else {
        wi <- t(xi) %*% mui_dev %*% vui %*% m0 %*% vui
        zi <- t(xi) %*% mui_dev %*% vui %*% m1 %*% vui
        
        ui0 <- (1/nsub)*wi %*% (yi-mui)
        ui1 <- (1/nsub)*zi %*% (yi-mui)
        
        ui[1:np,] <- ui0
        ui[(np+1):(2*np),] <- ui1
        
        C <- C + ui %*% t(ui)
        U <- U + ui
        
        if (is.na(C[1,1])) {
          print(iteration) # Commented out printing outside of print() methods
          print(ui)
          print(C)
        }
        
        di0 <- -(1/nsub) * wi %*% mui_dev %*% xi
        di1 <- -(1/nsub) * zi %*% mui_dev %*% xi
        
        firstdev[1:np,] <- di0
        firstdev[(np+1):(2*np),] <- di1
        Ufirstdev <- Ufirstdev + firstdev
      }
    }
    
    # calculate Q, betanew,
    Cinv=ginv(C)
    
    Q <- t(U) %*% Cinv %*% U
    
    arqif1dev <- t(Ufirstdev) %*% Cinv %*% U
    arqif2dev <- t(Ufirstdev) %*% Cinv %*% Ufirstdev
    
    invarqif2dev <- ginv(arqif2dev)
    
    betanew <- beta - invarqif2dev %*% arqif1dev
    betadiff <- abs(sum(betanew - beta))
    iteration <- iteration +1
  }
  
  
  ################  GMM end   ###########
  
  ################ Output 
  fit <- list()
  fit$terms <- Terms
  fit$formula <- as.vector(attr(Terms, "formula"))
  fit$call <- call
  fit$nobs <- nobs
  fit$iteration <- iteration
  fit$coefficients <- as.vector(beta)
  names(fit$coefficients) <- xnames
  fit$family <- family
  fit$y <- y
  fit$x <- x
  fit$id <- unique(id)
  fit$max.id <- max(nobs)
  fit$xnames <- xnames
  
  # covariance matrix
  dimnames(invarqif2dev)[[1]] <- xnames
  dimnames(invarqif2dev)[[2]] <- xnames
  fit$covariance <- invarqif2dev
  fit
}


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
b$event<- as.factor(b$event)
  
### Fit a GMM model ###

fit_GMM <- gmm_cloglog(pseudo ~ as.factor(Zb)+ as.factor(tpseudo),
                       id = id, data = b, family=gaussian, corstr="independence")
  
as.data.frame(cbind(mean = fit_GMM$coefficients, 
                    SE = sqrt(diag(fit_GMM$covariance))))
