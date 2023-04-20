## Nonlinear PK mixed-effects model, curved exponential family 
## Properties of the stochastic approximation algorithm for the evaluation of the two estimators of the Fisher information matrix
## (Section 4.2.1)

rm(list=ls())

## SAEM function without reprojection

saem <- function(data, nbiterem, nbiterburnin, theta0, kRW=0.5) {
  
  # data         : dataset 
  # nbiterem     : total number of iterations of the saem algorithm
  # nbiterburnin : number of burn-in iterations of the algorithm
  # theta0       : initial parameter values
  # kRW          : coefficient used to adjust the variance of the proposal kernel of the MCMC procedure
  
  # data processing
  xidep <- cbind(data$dose,data$time)
  y     <- data$y
  id    <- as.matrix(data$subject)
  n     <- length(unique(id))
  j     <- length(unique(data$time))
  
  
  # Model function
  
  model1cpt<-function(psi,id,xidep) { 
    dose  <- xidep[,1]
    tim   <- xidep[,2]  
    ka    <- psi[id,1]
    V     <- psi[id,2]
    CL    <- psi[id,3]
    k     <- CL/V
    ypred <- dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
    return(ypred)
  }
  
  # initial parameter values
  
  vpop     <- theta0$vpop
  kapop    <- theta0$kapop
  clpop    <- theta0$clpop
  omega2v  <- theta0$omega2v
  omega2ka <- theta0$omega2ka
  omega2cl <- theta0$omega2cl
  sigma2   <- theta0$sigma2
  
  p <- length(theta0)
  
  thetaest     <-  matrix(0,p,nbiterem)
  thetaest[,1] <- c(kapop,vpop,clpop,omega2ka,omega2v,omega2cl,sigma2)
  
  
  # variances of the proposal kernels of the MCMC procedure
  eta2v  <- kRW*omega2v
  eta2ka <- kRW*omega2ka
  eta2cl <- kRW*omega2cl
  
  # sequence of step sizes
  gamma <-c(rep(0.95,nbiterburnin), 1/(2:nbiterem)^0.6)
  
  # intermediary R objects
  deltaindi      <- array(0,c(p,n,nbiterem))
  H              <- array(0,c(p,p,n,nbiterem))
  G2             <- array(0,c(p,p,nbiterem))
  tempderiveeas  <-matrix(0,p,n)
  tempderiveeas2 <-matrix(0,p,p)
  dimstatexh     <- 7
  statexh        <- matrix(0,dimstatexh,nbiterem)

  
  # initial values for the individual parameters
  
  currentv   <- log(vpop)  + rnorm(n,0,sqrt(eta2v))
  currentka  <- log(kapop) + rnorm(n,0,sqrt(eta2ka))
  currentcl  <- log(clpop) + rnorm(n,0,sqrt(eta2cl))
  currentpsi <- cbind(exp(currentka),exp(currentv),exp(currentcl))
  
  ## Start of the em loop
  for (l in 1:(nbiterem-1)){
    
    ## Simulation step
    for (k in 1:(n)){
      
      ## Variable ka
      candidatka    <- currentka
      candidatka[k] <- candidatka[k] + rnorm(1,0,sqrt(eta2ka))
      psicandidat   <- cbind(exp(candidatka),exp(currentv),exp(currentcl))
      logs          <- -1/2/sigma2*sum((y-model1cpt(psicandidat,id,xidep))^2)+1/2/sigma2*sum((y-model1cpt(currentpsi,id,xidep))^2)
      logs          <- logs-1/2/omega2ka*((candidatka[k]-log(kapop))^2-(currentka[k]-log(kapop))^2)
      u             <- runif(1)
      logu          <- log(u)
      ind           <- (logu<logs)
      currentpsi    <- psicandidat*ind+currentpsi*(1-ind)
      currentka     <- candidatka*ind+currentka*(1-ind)
      
      ## Variable V
      candidatv    <- currentv
      candidatv[k] <- candidatv[k] + rnorm(1,0,sqrt(eta2v))
      psicandidat  <- cbind(exp(currentka),exp(candidatv),exp(currentcl))
      logs         <- -1/2/sigma2*sum((y-model1cpt(psicandidat,id,xidep))^2)+ 1/2/sigma2*sum((y-model1cpt(currentpsi,id,xidep))^2)
      logs         <- logs -1/2/omega2v*((candidatv[k]-log(vpop))^2-(currentv[k]-log(vpop))^2)
      u            <- runif(1)
      logu         <- log(u)
      ind          <- (logu<logs)
      currentpsi   <- psicandidat*ind+currentpsi*(1-ind)
      currentv     <- candidatv*ind+currentv*(1-ind)
      
      ## Variable cl
      candidatcl    <- currentcl
      candidatcl[k] <- candidatcl[k] + rnorm(1,0,sqrt(eta2cl))
      psicandidat   <- cbind(exp(currentka),exp(currentv),exp(candidatcl))
      logs          <- -1/2/sigma2*sum((y-model1cpt(psicandidat,id,xidep))^2)+ 1/2/sigma2*sum((y-model1cpt(currentpsi,id,xidep))^2)
      logs          <- logs -1/2/omega2cl*((candidatcl[k]-log(clpop))^2-(currentcl[k]-log(clpop))^2)
      u             <- runif(1)
      logu          <- log(u)
      ind           <- (logu<logs)
      currentpsi    <- psicandidat*ind+currentpsi*(1-ind)
      currentcl     <- candidatcl*ind+currentcl*(1-ind)
      
    }
    
    psi <- currentpsi
    
    # stochastic approximation of exhaustive statistics and parameter estimation update
    
    mco           <- matrix((y-model1cpt(psi,id,xidep))^2,n,j,byrow=TRUE)
    mcos          <- apply(mco,1,sum)
    STATEXH       <- c(apply(log(psi),2,mean), apply(log(psi)^2,2,mean), sum(mcos))
    statexh[,l+1] <- statexh[,l]*(1-gamma[l])+gamma[l]*STATEXH
    
       
    kapop    <- exp(statexh[1,l+1])
    vpop     <- exp(statexh[2,l+1])
    clpop    <- exp(statexh[3,l+1])
    omega2ka <- statexh[4,l+1]-statexh[1,l+1]^2
    omega2v  <- statexh[5,l+1]-statexh[2,l+1]^2
    omega2cl <- statexh[6,l+1]-statexh[3,l+1]^2
    sigma2   <- statexh[7,l+1]/n/j
    
    thetaest[,l+1] <- c(kapop, vpop, clpop, omega2ka, omega2v, omega2cl, sigma2)
    
    eta2ka   <- kRW*omega2ka
    eta2v    <- kRW*omega2v
    eta2cl   <- kRW*omega2cl

    ## Stochastic approximation of the derivatives of the complete log-likelihood
    
 
    tempderiveeas[1,]<- (log(psi[,1])-log(kapop))/omega2ka/kapop
    tempderiveeas[2,]<- (log(psi[,2])-log(vpop))/omega2v/vpop
    tempderiveeas[3,]<- (log(psi[,3])-log(clpop))/omega2cl/clpop
    tempderiveeas[4,]<- -1/2/omega2ka +1/2/omega2ka^2*(log(psi[,1])-log(kapop))^2
    tempderiveeas[5,]<- -1/2/omega2v +1/2/omega2v^2*(log(psi[,2])-log(vpop))^2
    tempderiveeas[6,]<- -1/2/omega2cl +1/2/omega2cl^2*(log(psi[,3])-log(clpop))^2
    tempderiveeas[7,]<- -j/2/sigma2+apply(mco,1,sum)/2/sigma2^2
    
    
    tempderiveeas2[1,1] <-  sum(-log(psi[,1])+log(kapop)-1)/omega2ka/kapop^2
    tempderiveeas2[2,2] <-  sum(-log(psi[,2])+log(vpop)-1)/omega2v/vpop^2
    tempderiveeas2[3,3] <-  sum(-log(psi[,3])+log(clpop)-1)/omega2cl/clpop^2
    tempderiveeas2[4,4] <-  n/2/omega2ka^2 -1/omega2ka^3*sum((log(psi[,1])-log(kapop))^2)
    tempderiveeas2[5,5] <-  n/2/omega2v^2 -1/omega2v^3*sum((log(psi[,2])-log(vpop))^2)
    tempderiveeas2[6,6] <-  n/2/omega2cl^2 -1/omega2cl^3*sum((log(psi[,3])-log(clpop))^2)
    tempderiveeas2[7,7] <-  n*j/2/sigma2^2-sum((y-model1cpt(psi,id,xidep))^2)/sigma2^3
    tempderiveeas2[1,4] <- -sum(log(psi[,1])-log(kapop))/omega2ka^2/kapop
    tempderiveeas2[2,5] <- -sum(log(psi[,2])-log(vpop))/omega2v^2/vpop
    tempderiveeas2[3,6] <- -sum(log(psi[,3])-log(clpop))/omega2cl^2/clpop
    tempderiveeas2[4,1] <- tempderiveeas2[1,4]
    tempderiveeas2[5,2] <- tempderiveeas2[2,5]
    tempderiveeas2[6,3] <- tempderiveeas2[3,6]
    
    
    deltaindi[,,l+1] <- deltaindi[,,l]*(1-gamma[l])+gamma[l]*tempderiveeas
    for (i in 1:n){
      H[,,i,l+1]<-H[,,i,l]*(1-gamma[l])+gamma[l]*(tempderiveeas[,i]%*%t(tempderiveeas[,i]))
    }
    G2[,,l+1] <- G2[,,l]*(1-gamma[l])+gamma[l]*(tempderiveeas2/n)
  }
  ## End of the em loop
  
  ## Computation of the FIM estimations
  
  isco <- array(0,c(p,p,nbiterem)) # estimation based on the score
  iobs <- array(0,c(p,p,nbiterem)) # estimation based on the observed information matrix
  SH   <- vector("list",nbiterem)
  
  
  for (t in 1:nbiterem){
    isco[,,t] <- deltaindi[,,t]%*%t(deltaindi[,,t])/n
    SH[[t]]<-matrix(0,p,p)
    for (i in 1:n){
      SH[[t]]<-SH[[t]]+H[,,i,t]
    }
    iobs[,,t] <- -G2[,,t] - SH[[t]]/n + isco[,,t]
  }
  
  res <- list(thetaest = thetaest, isco = isco, iobs = iobs)

  return(res)
}

## Definition of the sequences of compact sets

compact <- function(s,kappa){
  
  res <- (s[1]<=(20+kappa))*(s[1]>=(-20-kappa))*
    (s[2]<=(20+kappa))*(s[2]>=(-20-kappa))*
    (s[3]<=(20+kappa))*(s[3]>=(-20-kappa))*
    (s[4]<=(20+kappa))*
    (s[5]<=(20+kappa))*
    (s[6]<=(20+kappa))*
    (s[7]<=(50000))
  res <- as.numeric(res)
  return(res)
}

## SAEM function with reprojection

saem.proj <- function(data, nbiterem, nbiterburnin, theta0, kRW=0.5) {
  
  # data         : dataset 
  # nbiterem     : total number of iterations of the saem algorithm
  # nbiterburnin : number of burn-in iterations of the algorithm
  # theta0       : initial parameter values
  # kRW          : coefficient used to adjust the variance of the proposal kernel of the MCMC procedure
  
  
  
  # data processing
  xidep <- cbind(data$dose,data$time)
  y     <- data$y
  id    <- as.matrix(data$subject)
  n     <- length(unique(id))
  j     <- length(unique(data$time))
  
  
  # Model function and informations about the model
  
  model1cpt<-function(psi,id,xidep) { 
    dose  <- xidep[,1]
    tim   <- xidep[,2]  
    ka    <- psi[id,1]
    V     <- psi[id,2]
    CL    <- psi[id,3]
    k     <- CL/V
    ypred <- dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
    return(ypred)
  }
  
  nb.psi         <- 3
  dimstatexh     <- 7
  
  # initial parameter values
  
  vpop     <- theta0$vpop
  kapop    <- theta0$kapop
  clpop    <- theta0$clpop
  omega2v  <- theta0$omega2v
  omega2ka <- theta0$omega2ka
  omega2cl <- theta0$omega2cl
  sigma2   <- theta0$sigma2
  
  p <- length(theta0)
  
  thetaest     <-  matrix(0,p,nbiterem)
  thetaest[,1] <- c(kapop,vpop,clpop,omega2ka,omega2v,omega2cl,sigma2)
  
  
  # variances of the proposal kernels of the MCMC procedure
  eta2v  <- kRW*omega2v
  eta2ka <- kRW*omega2ka
  eta2cl <- kRW*omega2cl
  
  # sequence of step sizes
  gamma <-c(rep(0.95,nbiterburnin), 1/(2:nbiterem)^0.6)
  
  C       <- 50000
  eta     <- 2/5
  epsilon <- C*gamma^eta
  
  # initialize counters for reprojections
  
  nu    <- 1
  zeta  <- 1
  kappa <- 1
  
  # intermediary R objects
  deltaindi      <- array(0,c(p,n,nbiterem))
  H              <- array(0,c(p,p,n,nbiterem))
  G2             <- array(0,c(p,p,nbiterem))
  tempderiveeas  <- matrix(0,p,n)
  tempderiveeas2 <- matrix(0,p,p)
  statexh        <- matrix(0,dimstatexh,nbiterem)
  psi            <- array(0,dim=c(n,nb.psi,nbiterem))
  
  
  # initial values for the individual parameters and exhaustive statistics
  
  currentv   <- log(vpop)  + rnorm(n,0,sqrt(eta2v))
  currentka  <- log(kapop) + rnorm(n,0,sqrt(eta2ka))
  currentcl  <- log(clpop) + rnorm(n,0,sqrt(eta2cl))
  currentpsi <- cbind(exp(currentka),exp(currentv),exp(currentcl))
  psi[,,1]   <- psiinit <- currentpsi 
  mco           <- matrix((y-model1cpt(psi[,,1],id,xidep))^2,n,j,byrow=TRUE)
  mcos          <- apply(mco,1,sum)
  statexh[,1]       <- c(apply(log(psi[,,1]),2,mean), apply(log(psi[,,1])^2,2,mean), sum(mcos))
  
  while(compact(statexh[,1],kappa)==0){
    currentv   <- log(vpop)  + rnorm(n,0,sqrt(eta2v))
    currentka  <- log(kapop) + rnorm(n,0,sqrt(eta2ka))
    currentcl  <- log(clpop) + rnorm(n,0,sqrt(eta2cl))
    currentpsi <- cbind(exp(currentka),exp(currentv),exp(currentcl))
    psi[,,1]   <- psiinit <- currentpsi 
    mco        <- matrix((y-model1cpt(psi[,,1],id,xidep))^2,n,j,byrow=TRUE)
    mcos       <- apply(mco,1,sum)
    statexh[,1]<- c(apply(log(psi[,,1]),2,mean), apply(log(psi[,,1])^2,2,mean), sum(mcos))
  }
  
  ## Start of the em loop
  for (l in 1:(nbiterem-1)){
    
    
    ## Simulation step
    
    if (nu==0){ 
      ## Reprojection
      psi[,,l+1] <- currentpsi <- psiinit
      
    } else{
      ## Standard simulation procedure 
      for (k in 1:(n)){
        
        ## Variable ka
        candidatka    <- currentka
        candidatka[k] <- candidatka[k] + rnorm(1,0,sqrt(eta2ka))
        psicandidat   <- cbind(exp(candidatka),exp(currentv),exp(currentcl))
        logs          <- -1/2/sigma2*sum((y-model1cpt(psicandidat,id,xidep))^2)+1/2/sigma2*sum((y-model1cpt(currentpsi,id,xidep))^2)
        logs          <- logs-1/2/omega2ka*((candidatka[k]-log(kapop))^2-(currentka[k]-log(kapop))^2)
        u             <- runif(1)
        logu          <- log(u)
        ind           <- (logu<logs)
        currentpsi    <- psicandidat*ind+currentpsi*(1-ind)
        currentka     <- candidatka*ind+currentka*(1-ind)
        
        ## Variable V
        candidatv    <- currentv
        candidatv[k] <- candidatv[k] + rnorm(1,0,sqrt(eta2v))
        psicandidat  <- cbind(exp(currentka),exp(candidatv),exp(currentcl))
        logs         <- -1/2/sigma2*sum((y-model1cpt(psicandidat,id,xidep))^2)+ 1/2/sigma2*sum((y-model1cpt(currentpsi,id,xidep))^2)
        logs         <- logs -1/2/omega2v*((candidatv[k]-log(vpop))^2-(currentv[k]-log(vpop))^2)
        u            <- runif(1)
        logu         <- log(u)
        ind          <- (logu<logs)
        currentpsi   <- psicandidat*ind+currentpsi*(1-ind)
        currentv     <- candidatv*ind+currentv*(1-ind)
        
        ## Variable cl
        candidatcl    <- currentcl
        candidatcl[k] <- candidatcl[k] + rnorm(1,0,sqrt(eta2cl))
        psicandidat   <- cbind(exp(currentka),exp(currentv),exp(candidatcl))
        logs          <- -1/2/sigma2*sum((y-model1cpt(psicandidat,id,xidep))^2)+ 1/2/sigma2*sum((y-model1cpt(currentpsi,id,xidep))^2)
        logs          <- logs -1/2/omega2cl*((candidatcl[k]-log(clpop))^2-(currentcl[k]-log(clpop))^2)
        u             <- runif(1)
        logu          <- log(u)
        ind           <- (logu<logs)
        currentpsi    <- psicandidat*ind+currentpsi*(1-ind)
        currentcl     <- candidatcl*ind+currentcl*(1-ind)
        
      }
      psi[,,l+1] <- currentpsi
    }
    
    # stochastic approximation of exhaustive statistics and parameter estimation update
    
    mco           <- matrix((y-model1cpt(psi[,,l+1],id,xidep))^2,n,j,byrow=TRUE)
    mcos          <- apply(mco,1,sum)
    STATEXH       <- c(apply(log(psi[,,l+1]),2,mean), apply(log(psi[,,l+1])^2,2,mean), sum(mcos))
    statexh[,l+1] <- statexh[,l]*(1-gamma[l])+gamma[l]*STATEXH
    
    norm.delta.statexh <- sqrt(sum((statexh[,l+1]-statexh[,l])^2))
    
    ## check if reprojection will be necessary at next iteration
    if ((norm.delta.statexh<=epsilon[zeta]) && (compact(statexh[,l+1],kappa)==1)){
      kappa <- kappa
      zeta  <- zeta + 1
      nu    <- nu + 1
    } else{
      nu    <- 0
      kappa <- kappa + 1
      zeta  <- zeta + 1 
    }
    
    kapop    <- exp(statexh[1,l+1])
    vpop     <- exp(statexh[2,l+1])
    clpop    <- exp(statexh[3,l+1])
    omega2ka <- statexh[4,l+1]-statexh[1,l+1]^2
    omega2v  <- statexh[5,l+1]-statexh[2,l+1]^2
    omega2cl <- statexh[6,l+1]-statexh[3,l+1]^2
    sigma2   <- statexh[7,l+1]/n/j
    
    thetaest[,l+1] <- c(kapop, vpop, clpop, omega2ka, omega2v, omega2cl, sigma2)
    
    eta2ka   <- kRW*omega2ka
    eta2v    <- kRW*omega2v
    eta2cl   <- kRW*omega2cl
    
    ## Stochastic approximation of the derivatives of the complete log-likelihood
    
    
    tempderiveeas[1,]<- (log(psi[,1,l+1])-log(kapop))/omega2ka/kapop
    tempderiveeas[2,]<- (log(psi[,2,l+1])-log(vpop))/omega2v/vpop
    tempderiveeas[3,]<- (log(psi[,3,l+1])-log(clpop))/omega2cl/clpop
    tempderiveeas[4,]<- -1/2/omega2ka +1/2/omega2ka^2*(log(psi[,1,l+1])-log(kapop))^2
    tempderiveeas[5,]<- -1/2/omega2v +1/2/omega2v^2*(log(psi[,2,l+1])-log(vpop))^2
    tempderiveeas[6,]<- -1/2/omega2cl +1/2/omega2cl^2*(log(psi[,3,l+1])-log(clpop))^2
    tempderiveeas[7,]<- -j/2/sigma2+apply(mco,1,sum)/2/sigma2^2
    
    
    tempderiveeas2[1,1] <-  sum(-log(psi[,1,l+1])+log(kapop)-1)/omega2ka/kapop^2
    tempderiveeas2[2,2] <-  sum(-log(psi[,2,l+1])+log(vpop)-1)/omega2v/vpop^2
    tempderiveeas2[3,3] <-  sum(-log(psi[,3,l+1])+log(clpop)-1)/omega2cl/clpop^2
    tempderiveeas2[4,4] <-  n/2/omega2ka^2 -1/omega2ka^3*sum((log(psi[,1,l+1])-log(kapop))^2)
    tempderiveeas2[5,5] <-  n/2/omega2v^2 -1/omega2v^3*sum((log(psi[,2,l+1])-log(vpop))^2)
    tempderiveeas2[6,6] <-  n/2/omega2cl^2 -1/omega2cl^3*sum((log(psi[,3,l+1])-log(clpop))^2)
    tempderiveeas2[7,7] <-  n*j/2/sigma2^2-sum((y-model1cpt(psi[,,l+1],id,xidep))^2)/sigma2^3
    tempderiveeas2[1,4] <- -sum(log(psi[,1,l+1])-log(kapop))/omega2ka^2/kapop
    tempderiveeas2[2,5] <- -sum(log(psi[,2,l+1])-log(vpop))/omega2v^2/vpop
    tempderiveeas2[3,6] <- -sum(log(psi[,3,l+1])-log(clpop))/omega2cl^2/clpop
    tempderiveeas2[4,1] <- tempderiveeas2[1,4]
    tempderiveeas2[5,2] <- tempderiveeas2[2,5]
    tempderiveeas2[6,3] <- tempderiveeas2[3,6]
    
    
    deltaindi[,,l+1] <- deltaindi[,,l]*(1-gamma[l])+gamma[l]*tempderiveeas
    for (i in 1:n){
      H[,,i,l+1]<-H[,,i,l]*(1-gamma[l])+gamma[l]*(tempderiveeas[,i]%*%t(tempderiveeas[,i]))
    }
    G2[,,l+1] <- G2[,,l]*(1-gamma[l])+gamma[l]*(tempderiveeas2/n)
  }
  ## End of the em loop
  
  ## Computation of the FIM estimations
  
  isco <- array(0,c(p,p,nbiterem)) # estimation based on the score
  iobs <- array(0,c(p,p,nbiterem)) # estimation based on the observed information matrix
  SH   <- vector("list",nbiterem)
  
  
  for (t in 1:nbiterem){
    isco[,,t] <- deltaindi[,,t]%*%t(deltaindi[,,t])/n
    SH[[t]]<-matrix(0,p,p)
    for (i in 1:n){
      SH[[t]]<-SH[[t]]+H[,,i,t]
    }
    iobs[,,t] <- -G2[,,t] - SH[[t]]/n + isco[,,t]
  }
  
  res <- list(thetaest = thetaest, isco = isco, iobs = iobs)
  
  return(res)
}

## Monte-carlo estimation of the FIM

FIM_mc <- function(data, nbMC, nbMCburnin, theta, kRW=0.5) {
  
  # data       : dataset 
  # nbMC       : total number of iterations of Monte-Carlo iterations
  # nbMCburnin : number of burn-in iterations 
  # theta      : parameter values
  # kRW        : coefficient used to adjust the variance of the proposal kernel of the MCMC procedure
  
  
  # data processing
  xidep <- cbind(data$dose,data$time)
  y     <- data$y
  id    <- as.matrix(data$subject)
  n     <- length(unique(id))
  j     <- length(unique(data$time))
  
  # Model function
  
  model1cpt<-function(psi,id,xidep) { 
    dose  <- xidep[,1]
    tim   <- xidep[,2]  
    ka    <- psi[id,1]
    V     <- psi[id,2]
    CL    <- psi[id,3]
    k     <- CL/V
    ypred <- dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
    return(ypred)
  }
  
  # parameter values 
  
  vpop     <- theta$vpop
  kapop    <- theta$kapop
  clpop    <- theta$clpop
  omega2v  <- theta$omega2v
  omega2ka <- theta$omega2ka
  omega2cl <- theta$omega2cl
  sigma2   <- theta$sigma2
  
  p <- length(theta)
  
  # variances of the proposal kernels of the MCMC procedure
  eta2v  <- kRW*omega2v
  eta2ka <- kRW*omega2ka
  eta2cl <- kRW*omega2cl
  
  # intermediary R quantities
  psiMC          <- array(0,c(n,3,nbMC))
  iscoMC         <- array(0,c(p,p,nbMC))
  iobsMC         <- array(0,c(p,p,nbMC))
  temp2MC        <- array(0,c(p,p,nbMC))
  temp3MC        <- array(0,c(p,p,nbMC))
  tempderiveeMC  <- array(0,c(p,n,nbMC))
  tempderiveeMC2 <- array(0,c(p,p,nbMC))
  
  # initial values for the individual parameters
  currentv   <- log(vpop) + rnorm(n,0,sqrt(eta2v))
  currentka  <- log(kapop) + rnorm(n,0,sqrt(eta2ka))
  currentcl  <- log(clpop) + rnorm(n,0,sqrt(eta2cl))
  currentpsi <- cbind(exp(currentka),exp(currentv),exp(currentcl))
  
  # simulate a markov chain with stationnary distribution p(psi|y) 
  
  # burn-in
  
  for (m in 1:nbMCburnin){
    
    for (k in 1:n){
      
      # variable ka 
      
      candidatka    <- currentka
      candidatka[k] <- candidatka[k] + rnorm(1,0,sqrt(eta2ka))
      psicandidat   <- cbind(exp(candidatka),exp(currentv),exp(currentcl))
      logs          <- -1/2/sigma2*sum((y-model1cpt(psicandidat,id,xidep))^2)+1/2/sigma2*sum((y-model1cpt(currentpsi,id,xidep))^2)
      logs          <- logs-1/2/omega2ka*((candidatka[k]-log(kapop))^2-(currentka[k]-log(kapop))^2)
      u             <- runif(1)
      logu          <- log(u)
      ind           <- (logu<logs)
      currentpsi    <- psicandidat*ind+currentpsi*(1-ind)
      currentka     <- candidatka*ind+currentka*(1-ind)
      
      # variable v 
      
      candidatv     <- currentv
      candidatv[k]  <- candidatv[k] + rnorm(1,0,sqrt(eta2v))
      psicandidat   <- cbind(exp(currentka),exp(candidatv),exp(currentcl))
      logs          <- -1/2/sigma2*sum((y-model1cpt(psicandidat,id,xidep))^2)+1/2/sigma2*sum((y-model1cpt(currentpsi,id,xidep))^2)
      logs          <- logs-1/2/omega2v*((candidatv[k]-log(vpop))^2-(currentv[k]-log(vpop))^2)
      u             <- runif(1)
      logu          <- log(u)
      ind           <- (logu<logs)
      currentpsi    <- psicandidat*ind+currentpsi*(1-ind)
      currentv      <- candidatv*ind+currentv*(1-ind)
      
      # variable cl 
      
      candidatcl    <- currentcl
      candidatcl[k] <- candidatcl[k] + rnorm(1,0,sqrt(eta2cl))
      psicandidat   <- cbind(exp(currentka),exp(currentv),exp(candidatcl))
      logs          <- -1/2/sigma2*sum((y-model1cpt(psicandidat,id,xidep))^2)+1/2/sigma2*sum((y-model1cpt(currentpsi,id,xidep))^2)
      logs          <- logs-1/2/omega2cl*((candidatcl[k]-log(clpop))^2-(currentcl[k]-log(clpop))^2)
      u             <- runif(1)
      logu          <- log(u)
      ind           <- (logu<logs)
      currentpsi    <- psicandidat*ind+currentpsi*(1-ind)
      currentcl     <- candidatcl*ind+currentcl*(1-ind)
      
    }
  }
  
  
  for (m in 1:nbMC){
    
    for (k in 1:n){
      
      # variable ka 
      
      candidatka    <- currentka
      candidatka[k] <- candidatka[k] + rnorm(1,0,sqrt(eta2ka))
      psicandidat   <- cbind(exp(candidatka),exp(currentv),exp(currentcl))
      logs          <- -1/2/sigma2*sum((y-model1cpt(psicandidat,id,xidep))^2)+1/2/sigma2*sum((y-model1cpt(currentpsi,id,xidep))^2)
      logs          <- logs-1/2/omega2ka*((candidatka[k]-log(kapop))^2-(currentka[k]-log(kapop))^2)
      u             <- runif(1)
      logu          <- log(u)
      ind           <- (logu<logs)
      currentpsi    <- psicandidat*ind+currentpsi*(1-ind)
      currentka     <- candidatka*ind+currentka*(1-ind)
      
      # variable v 
      
      candidatv     <- currentv
      candidatv[k]  <-candidatv[k] + rnorm(1,0,sqrt(eta2v))
      psicandidat   <- cbind(exp(currentka),exp(candidatv),exp(currentcl))
      logs          <- -1/2/sigma2*sum((y-model1cpt(psicandidat,id,xidep))^2)+1/2/sigma2*sum((y-model1cpt(currentpsi,id,xidep))^2)
      logs          <- logs-1/2/omega2v*((candidatv[k]-log(vpop))^2-(currentv[k]-log(vpop))^2)
      u             <- runif(1)
      logu          <- log(u)
      ind           <- (logu<logs)
      currentpsi    <- psicandidat*ind+currentpsi*(1-ind)
      currentv      <- candidatv*ind+currentv*(1-ind)
      
      # variable cl 
      
      candidatcl    <- currentcl
      candidatcl[k] <- candidatcl[k] + rnorm(1,0,sqrt(eta2cl))
      psicandidat   <- cbind(exp(currentka),exp(currentv),exp(candidatcl))
      logs          <- -1/2/sigma2*sum((y-model1cpt(psicandidat,id,xidep))^2)+1/2/sigma2*sum((y-model1cpt(currentpsi,id,xidep))^2)
      logs          <- logs-1/2/omega2cl*((candidatcl[k]-log(clpop))^2-(currentcl[k]-log(clpop))^2)
      u             <- runif(1)
      logu          <- log(u)
      ind           <- (logu<logs)
      currentpsi    <- psicandidat*ind+currentpsi*(1-ind)
      currentcl     <- candidatcl*ind+currentcl*(1-ind)
      
    }
    psiMC[,,m] <- currentpsi
  }
  
  ## Monte-Carlo evaluation of the FIM estimator based on the score
  
  # first order derivatives of the individual complete log-likelihood
  tempderiveeMC[1,,] <- (log(psiMC[,1,])-log(kapop))/omega2ka/kapop
  tempderiveeMC[2,,] <- (log(psiMC[,2,])-log(vpop))/omega2v/vpop
  tempderiveeMC[3,,] <- (log(psiMC[,3,])-log(clpop))/omega2cl/clpop
  tempderiveeMC[4,,] <- -1/2/omega2ka+1/2/omega2ka^2*(log(psiMC[,1,])-log(kapop))^2
  tempderiveeMC[5,,] <- -1/2/omega2v+1/2/omega2v^2*(log(psiMC[,2,])-log(vpop))^2
  tempderiveeMC[6,,] <- -1/2/omega2cl+1/2/omega2cl^2*(log(psiMC[,3,])-log(clpop))^2
  for (m in 1:nbMC){
    mco <- matrix((y-model1cpt(psiMC[,,m],id,xidep))^2,n,j,byrow=TRUE)
    tempderiveeMC[7,,m]<- -j/2/sigma2+apply(mco,1,sum)/2/sigma2^2
  }
  
  # computation of the FIM estimator based on the score (Isco)
  for (m in 1:nbMC)
  {
    SnbMCindi <- apply(tempderiveeMC[,,1:m],c(1,2),mean)
    temp <- 0
    for (i in 1:n){
      temp <- SnbMCindi[,i]%*%t(SnbMCindi[,i])+temp
    }
    iscoMC[,,m]<-temp/n
  }
  
  ## Monte-Carlo evaluation of the observed FIM according to the Louis' formula
  
  # second order derivatives of the individual complete log-likelihood
  for (m in 1:nbMC){
    tempderiveeMC2[1,1,m] <- sum(-log(psiMC[,1,m])+log(kapop)-1)/omega2ka/kapop^2
    tempderiveeMC2[2,2,m] <- sum(-log(psiMC[,2,m])+log(vpop)-1)/omega2v/vpop^2
    tempderiveeMC2[3,3,m] <- sum(-log(psiMC[,3,m])+log(clpop)-1)/omega2cl/clpop^2
    tempderiveeMC2[4,4,m] <- n/2/omega2ka^2 -1/omega2ka^3*sum((log(psiMC[,1,m])-log(kapop))^2)
    tempderiveeMC2[5,5,m] <- n/2/omega2v^2 -1/omega2v^3*sum((log(psiMC[,2,m])-log(vpop))^2)
    tempderiveeMC2[6,6,m] <- n/2/omega2cl^2 -1/omega2cl^3*sum((log(psiMC[,3,m])-log(clpop))^2)
    tempderiveeMC2[7,7,m] <- n*j/2/sigma2^2-sum((y-model1cpt(psiMC[,,m],id,xidep))^2)/sigma2^3
    tempderiveeMC2[1,4,m]<- -sum(log(psiMC[,1,m])-log(kapop))/omega2ka^2/kapop
    tempderiveeMC2[2,5,m]<- -sum(log(psiMC[,2,m])-log(vpop))/omega2v^2/vpop
    tempderiveeMC2[3,6,m]<- -sum(log(psiMC[,3,m])-log(clpop))/omega2cl^2/clpop
    tempderiveeMC2[4,1,m]<- tempderiveeMC2[1,4,m]
    tempderiveeMC2[5,2,m]<- tempderiveeMC2[2,5,m]
    tempderiveeMC2[6,3,m]<- tempderiveeMC2[3,6,m]
  }
  
  
  
  # computation of the observed FIM (Iobs)
  
  for (m in 1:nbMC){
    temp2MC[,,m]<-apply(tempderiveeMC2[,,1:m],c(1,2),sum)/m/n
  }
  temp3 <- 0
  for (m in 1:nbMC){
    for (i in 1:n){
      temp3 <-   tempderiveeMC[,i,m]%*%t(tempderiveeMC[,i,m])+temp3
    }
    temp3MC[,,m] <- temp3/m/n
  }
  TnbMC<- (temp2MC + temp3MC)
  iobsMC <- iscoMC - TnbMC
  
  res <- list(iobsMC = iobsMC[,,nbMC], iscoMC =iscoMC[,,nbMC]) 
  return(res)
  
}

## 1- Data simulation

# Model function
## returns the conditional mean of the observations given the random parameters psi and
## structural covariates, in this case dose and time


model1cptsim<-function(psi,id,tim,dose) { 
  ka   <-psi[id,1]
  V    <-psi[id,2]
  CL   <-psi[id,3]
  k    <-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

# Sample characteristics  

n     <- 100                             # number of subjects
times <- c(0.25,0.5,1,2,3.5,5,7,9,12,24) # observation times
j     <- length(times)                   # number of observations per subject
dose  <- 320                             # dose

# True parameter values

vpop     <- 31
kapop    <- 1.6
clpop    <- 2.8
omega2v  <- 0.40
omega2ka <- 0.40
omega2cl <- 0.40
sigma2   <- 0.75

# Simulation of the individual parameters

vind  <- exp(rnorm(n,log(vpop),sd=sqrt(omega2v)))
kaind <- exp(rnorm(n,log(kapop),sd=sqrt(omega2ka)))
clind <- exp(rnorm(n,log(clpop),sd=sqrt(omega2cl)))

# Simulation of the observations

ypred <- c()

for (k in 1:n){
  ypred <- c(ypred,model1cptsim(cbind(kaind,vind,clind),k, times,dose))
}

y <- ypred + rnorm(n*j,0,sd=sqrt(sigma2))

datasim <- data.frame(y=y,dose=rep(dose,n*j),time=rep(times,n),subject=kronecker(1:n, rep(1,j)))


## 2- Numerical experiment on the same simulated dataset

## a- Evaluation of both estimators of the FIM using the stochastic approximation algorithm

nbsim <- 1 #number of replicates

# Algorithmic settings
nbiterem <- 3000
nbiterburnin <- 1000

# Saving the nbsim results
iscoarray <- array(0,dim=c(nbsim,7,7,nbiterem))
iobsarray <- array(0,dim=c(nbsim,7,7,nbiterem))
thetaest <- matrix(NA,7,nbsim)

for (k in 1:nbsim){
  
  set.seed(k*100+10)
  
  theta0         <- list(vpop=vpop*runif(1,0.8,1.2),kapop=kapop*runif(1,0.8,1.2),clpop=clpop*runif(1,0.8,1.2),omega2v=omega2v*runif(1,0.4,2),
                 omega2ka=omega2ka*runif(1,0.4,2),omega2cl=omega2cl*runif(1,0.4,2),sigma2=sigma2*runif(1,0.4,2))
  res            <- saem.proj(datasim, nbiterem, nbiterburnin, theta0)
  #res            <- saem(datasim, nbiterem, nbiterburnin, theta0)
  iscoarray[k,,,]<- res$isco 
  iobsarray[k,,,]<- res$iobs 
  thetaest[,k]   <- res$thetaest[,nbiterem]
}

## b- Monte-Carlo evaluation of both estimates, considered as the targets for the estimates computed using the stochastic approximation algorithm

nbMC <- 10000
nbMCburnin <- 5000

tm <- rowMeans(thetaest)

thetaMean <- list(kapop=tm[1],vpop=tm[2],clpop=tm[3],omega2ka=tm[4],omega2v=tm[5],omega2cl=tm[6],sigma2=tm[7])

FisherMC <- FIM_mc(datasim, nbMC, nbMCburnin, thetaMean)

iscoMC <- FisherMC$iscoMC
iobsMC <- FisherMC$iobsMC

###
### Generation of results graphs
###

library(ggplot2)
library(cowplot)
library(devtools)
library(RColorBrewer)

#i- Bias over iterations

biasIsco <- array(0,dim=c(7,7,nbsim,nbiterem))

for (j in 1:nbsim){
  for (k in 1:nbiterem){
    biasIsco[,,j,k] <- (iscoarray[j,,,k] - iscoMC)/iscoMC
  }
}

biasIobs <- array(0,dim=c(7,7,nbsim,nbiterem))

for (j in 1:nbsim){
  for (k in 1:nbiterem){
    biasIobs[,,j,k] <- (iobsarray[j,,,k] - iobsMC)/iobsMC
  }
}


DataResBias <- data.frame(
  BiasF11=c(abs(colMeans(biasIsco[1,1,,])[(nbiterburnin+1):nbiterem]),abs(colMeans(biasIobs[1,1,,])[(nbiterburnin+1):nbiterem])),
  BiasF22=c(abs(colMeans(biasIsco[2,2,,])[(nbiterburnin+1):nbiterem]),abs(colMeans(biasIobs[2,2,,])[(nbiterburnin+1):nbiterem])),
  BiasF33=c(abs(colMeans(biasIsco[3,3,,])[(nbiterburnin+1):nbiterem]),abs(colMeans(biasIobs[3,3,,])[(nbiterburnin+1):nbiterem])),
  BiasF44=c(abs(colMeans(biasIsco[4,4,,])[(nbiterburnin+1):nbiterem]),abs(colMeans(biasIobs[4,4,,])[(nbiterburnin+1):nbiterem])),
  BiasF55=c(abs(colMeans(biasIsco[5,5,,])[(nbiterburnin+1):nbiterem]),abs(colMeans(biasIobs[5,5,,])[(nbiterburnin+1):nbiterem])),
  BiasF66=c(abs(colMeans(biasIsco[6,6,,])[(nbiterburnin+1):nbiterem]),abs(colMeans(biasIobs[6,6,,])[(nbiterburnin+1):nbiterem])),
  Iter=c(seq(1,nbiterem-nbiterburnin),seq(1,nbiterem-nbiterburnin)),
  Estimate=c(rep('Isco',nbiterem-nbiterburnin),rep('Iobs',nbiterem-nbiterburnin)))

BiasF11 <- ggplot(DataResBias, aes(y=BiasF11, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~ka~','~ ka~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

BiasF22 <- ggplot(DataResBias, aes(y=BiasF22, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~V~','~ V~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

BiasF33 <- ggplot(DataResBias, aes(y=BiasF33, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~Cl~','~ Cl~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

BiasF44 <- ggplot(DataResBias, aes(y=BiasF44, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~omega[ka]^2~','~ omega[ka]^2~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

BiasF55 <- ggplot(DataResBias, aes(y=BiasF55, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~omega[V]^2~','~ omega[V]^2~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

BiasF66 <- ggplot(DataResBias, aes(y=BiasF66, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~omega[Cl]^2~','~ omega[Cl]^2~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))


plot_grid(BiasF11, BiasF22, BiasF33, BiasF44, BiasF55, BiasF66, ncol = 3, nrow = 2)

#ii- Relative standard deviations over iterations

rsdIsco <- array(0,dim=c(7,7,nbsim,nbiterem))

for (j in 1:nbsim){
  for (k in 1:nbiterem){
    rsdIsco[,,j,k] <- (iscoarray[j,,,k] - iscoMC)^2/iscoMC^2
  }
}

rsdIobs <- array(0,dim=c(7,7,nbsim,nbiterem))

for (j in 1:nbsim){
  for (k in 1:nbiterem){
    rsdIobs[,,j,k] <- (iobsarray[j,,,k] - iobsMC)^2/iobsMC^2
  }
}

# Graphes des relative standard deviations au fil des itérations

DataResRsd <- data.frame(
  RsdF11=c(colMeans(rsdIsco[1,1,,])[(nbiterburnin+1):nbiterem],colMeans(rsdIobs[1,1,,])[(nbiterburnin+1):nbiterem]),
  RsdF22=c(colMeans(rsdIsco[2,2,,])[(nbiterburnin+1):nbiterem],colMeans(rsdIobs[2,2,,])[(nbiterburnin+1):nbiterem]),
  RsdF33=c(colMeans(rsdIsco[3,3,,])[(nbiterburnin+1):nbiterem],colMeans(rsdIobs[3,3,,])[(nbiterburnin+1):nbiterem]),
  RsdF44=c(colMeans(rsdIsco[4,4,,])[(nbiterburnin+1):nbiterem],colMeans(rsdIobs[4,4,,])[(nbiterburnin+1):nbiterem]),
  RsdF55=c(colMeans(rsdIsco[5,5,,])[(nbiterburnin+1):nbiterem],colMeans(rsdIobs[5,5,,])[(nbiterburnin+1):nbiterem]),
  RsdF66=c(colMeans(rsdIsco[6,6,,])[(nbiterburnin+1):nbiterem],colMeans(rsdIobs[6,6,,])[(nbiterburnin+1):nbiterem]),
  Iter=c(seq(1,nbiterem-nbiterburnin),seq(1,nbiterem-nbiterburnin)),
  Estimate=c(rep('Vn',nbiterem-nbiterburnin),rep('Wn',nbiterem-nbiterburnin)))

RsdF11 <- ggplot(DataResRsd, aes(y=RsdF11, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~ka~','~ ka~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

RsdF22 <- ggplot(DataResRsd, aes(y=RsdF22, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~V~','~ V~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

RsdF33 <- ggplot(DataResRsd, aes(y=RsdF33, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~Cl~','~ Cl~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

RsdF44 <- ggplot(DataResRsd, aes(y=RsdF44, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~omega[ka]^2~','~ omega[ka]^2~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

RsdF55 <- ggplot(DataResRsd, aes(y=RsdF55, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~omega[V]^2~','~ omega[V]^2~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

RsdF66 <- ggplot(DataResRsd, aes(y=RsdF66, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~omega[Cl]^2~','~ omega[Cl]^2~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

plot_grid(RsdF11, RsdF22, RsdF33, RsdF44, RsdF55, RsdF66, ncol = 3, nrow = 2)
