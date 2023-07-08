## R function implementing the saem algorithm to compute the parameter estimates 
## and the FIM estimates simultaneously in the PK nonlinear mixed-effects model
## with three random effects, thus belonging to the curved exponential family

saem <- function(data, nbiterem, nbiterburnin, theta0, kRW=0.5) {
  
  # data         : dataset 
  # nbiterem     : total number of iterations of the saem algorithm
  # nbiterburnin : number of burn-in iterations of the algorithm
  # theta0       : initial parameter values
  # kRW          : coefficient used to adjust the variance of the proposal 
  #                kernel of the MCMC procedure
  
  # data processing
  xidep <- cbind(data$dose,data$time)
  y     <- data$y
  id    <- as.matrix(data$subject)
  n     <- length(unique(id))
  j     <- length(unique(data$time))
  
  
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
  tempderiveeas  <-matrix(0,p,n)
  tempderiveeas2 <-matrix(0,p,p)
  dimstatexh     <- 7
  statexh        <- array(0,dim=c(dimstatexh,n,nbiterem)) 
  psi            <- array(0,dim=c(n,nb.psi,nbiterem))
  
  # initial values for the individual parameters
  
  currentv   <- log(vpop)  + rnorm(n,0,sqrt(eta2v))
  currentka  <- log(kapop) + rnorm(n,0,sqrt(eta2ka))
  currentcl  <- log(clpop) + rnorm(n,0,sqrt(eta2cl))
  currentpsi <- cbind(exp(currentka),exp(currentv),exp(currentcl))
  
  psi[,,1]   <- psiinit <- currentpsi 
  
  while(compact(statexh[,,1],kappa)==0){
    currentv   <- log(vpop)  + rnorm(n,0,sqrt(eta2v))
    currentka  <- log(kapop) + rnorm(n,0,sqrt(eta2ka))
    currentcl  <- log(clpop) + rnorm(n,0,sqrt(eta2cl))
    currentpsi <- cbind(exp(currentka),exp(currentv),exp(currentcl))
    psi[,,1]   <- psiinit <- currentpsi 
    mco           <- matrix((y-model1cpt(psi,id,xidep))^2,n,j,byrow=TRUE)
    statexh[1:3,,1] <- t(log(psi))
    statexh[4:6,,1] <- t(log(psi)^2)
    statexh[7,,1]   <- apply(mco,1,sum)
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
        logs          <- -1/2/sigma2*sum((y-model1cpt(psicandidat,id,xidep))^2)+
          1/2/sigma2*sum((y-model1cpt(currentpsi,id,xidep))^2)
        logs          <- logs-1/2/omega2ka*((candidatka[k]-log(kapop))^2-
                                              (currentka[k]-log(kapop))^2)
        u             <- runif(1)
        logu          <- log(u)
        ind           <- (logu<logs)
        currentpsi    <- psicandidat*ind+currentpsi*(1-ind)
        currentka     <- candidatka*ind+currentka*(1-ind)
        
        ## Variable V
        candidatv    <- currentv
        candidatv[k] <- candidatv[k] + rnorm(1,0,sqrt(eta2v))
        psicandidat  <- cbind(exp(currentka),exp(candidatv),exp(currentcl))
        logs         <- -1/2/sigma2*sum((y-model1cpt(psicandidat,id,xidep))^2)+ 
          1/2/sigma2*sum((y-model1cpt(currentpsi,id,xidep))^2)
        logs         <- logs -
          1/2/omega2v*((candidatv[k]-log(vpop))^2-(currentv[k]-log(vpop))^2)
        u            <- runif(1)
        logu         <- log(u)
        ind          <- (logu<logs)
        currentpsi   <- psicandidat*ind+currentpsi*(1-ind)
        currentv     <- candidatv*ind+currentv*(1-ind)
        
        ## Variable cl
        candidatcl    <- currentcl
        candidatcl[k] <- candidatcl[k] + rnorm(1,0,sqrt(eta2cl))
        psicandidat   <- cbind(exp(currentka),exp(currentv),exp(candidatcl))
        logs          <- -1/2/sigma2*sum((y-model1cpt(psicandidat,id,xidep))^2)+ 
          1/2/sigma2*sum((y-model1cpt(currentpsi,id,xidep))^2)
        logs          <- logs -
          1/2/omega2cl*((candidatcl[k]-log(clpop))^2-(currentcl[k]-log(clpop))^2)
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
    statexh[1:3,,l+1] <- statexh[1:3,,l]*(1-gamma[l])+gamma[l]*t(log(psi[,,l+1]))
    statexh[4:6,,l+1] <- statexh[4:6,,l]*(1-gamma[l])+gamma[l]*t(log(psi[,,l+1])^2)
    statexh[7,,l+1]   <- statexh[7,,l]*(1-gamma[l])+gamma[l]*apply(mco,1,sum)
    
    
    norm.delta.statexh <- sqrt(sum((statexh[,,l+1]-statexh[,,l])^2))
    
    ## check if reprojection will be necessary at next iteration
    if ((norm.delta.statexh<=epsilon[zeta]) && (compact(statexh[,,l+1],kappa)==1)){
      kappa <- kappa
      zeta  <- zeta + 1
      nu    <- nu + 1
    } else{
      nu    <- 0
      kappa <- kappa + 1
      zeta  <- zeta + 1 
    }
    
    
    kapop    <- exp(mean(statexh[1,,l+1]))
    vpop     <- exp(mean(statexh[2,,l+1]))
    clpop    <- exp(mean(statexh[3,,l+1]))
    omega2ka <- mean(statexh[4,,l+1])-mean(statexh[1,,l+1])^2
    omega2v  <- mean(statexh[5,,l+1])-mean(statexh[2,,l+1])^2
    omega2cl <- mean(statexh[6,,l+1])-mean(statexh[3,,l+1])^2
    sigma2   <- sum(statexh[7,,l+1])/n/j
    
    thetaest[,l+1] <- c(kapop, vpop, clpop, omega2ka, omega2v, omega2cl, sigma2)
    
    eta2ka   <- kRW*omega2ka
    eta2v    <- kRW*omega2v
    eta2cl   <- kRW*omega2cl
    
    ## Stochastic approximation of the derivatives of the complete log-likelihood
    
    ### For the computation of Isco
    
    deltaindi[1,,l+1]<- (statexh[1,,l+1]-log(kapop))/omega2ka/kapop
    deltaindi[2,,l+1]<- (statexh[2,,l+1]-log(vpop))/omega2v/vpop
    deltaindi[3,,l+1]<- (statexh[3,,l+1]-log(clpop))/omega2cl/clpop
    deltaindi[4,,l+1]<- -1/2/omega2ka +
      1/2/omega2ka^2*(statexh[4,,l+1]-2*statexh[1,,l+1]*log(kapop)+log(kapop)^2)
    deltaindi[5,,l+1]<- -1/2/omega2v +
      1/2/omega2v^2*(statexh[5,,l+1]-2*statexh[2,,l+1]*log(vpop)+log(vpop)^2)
    deltaindi[6,,l+1]<- -1/2/omega2cl +
      1/2/omega2cl^2*(statexh[6,,l+1]-2*statexh[3,,l+1]*log(clpop)+log(clpop)^2)
    deltaindi[7,,l+1]<- -j/2/sigma2+statexh[7,,l+1]/2/sigma2^2
    
    ### For the computation of Iobs
    
    tempderiveeas[1,]<- (log(psi[,1,l+1])-log(kapop))/omega2ka/kapop
    tempderiveeas[2,]<- (log(psi[,2,l+1])-log(vpop))/omega2v/vpop
    tempderiveeas[3,]<- (log(psi[,3,l+1])-log(clpop))/omega2cl/clpop
    tempderiveeas[4,]<- -1/2/omega2ka +
      1/2/omega2ka^2*(log(psi[,1,l+1])-log(kapop))^2
    tempderiveeas[5,]<- -1/2/omega2v +
      1/2/omega2v^2*(log(psi[,2,l+1])-log(vpop))^2
    tempderiveeas[6,]<- -1/2/omega2cl +
      1/2/omega2cl^2*(log(psi[,3,l+1])-log(clpop))^2
    tempderiveeas[7,]<- -j/2/sigma2+apply(mco,1,sum)/2/sigma2^2
    
    
    tempderiveeas2[1,1] <-  sum(-log(psi[,1,l+1])+log(kapop)-1)/omega2ka/kapop^2
    tempderiveeas2[2,2] <-  sum(-log(psi[,2,l+1])+log(vpop)-1)/omega2v/vpop^2
    tempderiveeas2[3,3] <-  sum(-log(psi[,3,l+1])+log(clpop)-1)/omega2cl/clpop^2
    tempderiveeas2[4,4] <-  n/2/omega2ka^2 -
      1/omega2ka^3*sum((log(psi[,1,l+1])-log(kapop))^2)
    tempderiveeas2[5,5] <-  n/2/omega2v^2 -
      1/omega2v^3*sum((log(psi[,2,l+1])-log(vpop))^2)
    tempderiveeas2[6,6] <-  n/2/omega2cl^2 -
      1/omega2cl^3*sum((log(psi[,3,l+1])-log(clpop))^2)
    tempderiveeas2[7,7] <-  n*j/2/sigma2^2-
      sum((y-model1cpt(psi[,,l+1],id,xidep))^2)/sigma2^3
    tempderiveeas2[1,4] <- -sum(log(psi[,1,l+1])-log(kapop))/omega2ka^2/kapop
    tempderiveeas2[2,5] <- -sum(log(psi[,2,l+1])-log(vpop))/omega2v^2/vpop
    tempderiveeas2[3,6] <- -sum(log(psi[,3,l+1])-log(clpop))/omega2cl^2/clpop
    tempderiveeas2[4,1] <- tempderiveeas2[1,4]
    tempderiveeas2[5,2] <- tempderiveeas2[2,5]
    tempderiveeas2[6,3] <- tempderiveeas2[3,6]
    
    
    for (i in 1:n){
      H[,,i,l+1]<-H[,,i,l]*(1-gamma[l])+
        gamma[l]*(tempderiveeas[,i]%*%t(tempderiveeas[,i]))
    }
    G2[,,l+1] <- G2[,,l]*(1-gamma[l])+gamma[l]*(tempderiveeas2/n)
  }
  ## End of the em loop
  
  ## Computation of the FIM estimations
  
  isco <- array(0,c(p,p,nbiterem)) 
  iobs <- array(0,c(p,p,nbiterem)) 
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



compact <- function(s,kappa){
  
  res <- prod((s[1,]<=(20+kappa))*(s[1,]>=(-20-kappa))*
                (s[2,]<=(20+kappa))*(s[2,]>=(-20-kappa))*
                (s[3,]<=(20+kappa))*(s[3,]>=(-20-kappa))*
                (s[4,]<=(20+kappa))*
                (s[5,]<=(20+kappa))*
                (s[6,]<=(20+kappa))*
                (s[7,]<=(50000))) 
  res <- as.numeric(res)
  return(res)
}