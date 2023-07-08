## R function implementing the saem algorithm to compute the parameter estimates 
## and the FIM estimates simultaneously in the PK nonlinear mixed-effects model
## with two random effects, thus not belonging to the curved exponential family

saem_non_exp <- function(data, nbiterem, nbiterburnin, theta0, kRW=0.5) {
  
  # data         : dataset 
  # nbiterem     : total number of iterations of the SAEM algorithm
  # nbiterburnin : number of burn-in iterations of the algorithm
  # theta0       : initial parameter values
  # kRW          : coefficient used to adjust the variance of the proposal kernel 
  #                of the MCMC procedure
  
  # Q quantity 
  
  floglik <- function(v,y,psi,xidep,id,alpha){
    l <- length(alpha)
    psi <- array(psi,dim=c(n,2,l))
    value <- 0
    for (ll in 1:l){
      moyij <- model1cptV(psi[,,ll],id,xidep,v)
      value <- value + alpha[ll]*sum((y-moyij)^2) 
    }
    return(value)
  }
  
  # data processing
  xidep <- cbind(data$dose,data$time)
  y     <- data$y
  id    <- as.matrix(data$subject)
  n     <- length(unique(id))
  j     <- length(unique(data$time))
  
  # initial parameter values
  
  vpop     <- theta0$vpop
  kapop    <- theta0$kapop
  clpop    <- theta0$clpop
  omega2ka <- theta0$omega2ka
  omega2cl <- theta0$omega2cl
  sigma2   <- theta0$sigma2
  
  p <- length(theta0) 
  
  thetaest     <- matrix(0,p,nbiterem)
  thetaest[,1] <- c(kapop, vpop, clpop, omega2ka, omega2cl, sigma2)
  
  # variances of the proposal kernels of the MCMC procedure
  eta2ka  <- kRW*omega2ka
  eta2cl  <- kRW*omega2cl
  
  # sequence of step sizes 
  gamma <-  1/(1:(nbiterem))^0.501
  
  cumgamma <- matrix(0,nbiterem,nbiterem) 
  
  diag(cumgamma) <- gamma
  
  for (l in 2:nbiterem){
    for (m in 1:(l-1)){
      cumgamma[l,m] <- (1-gamma[l])*cumgamma[l-1,m]  
    }  
  }
  
  # intermediary R objects
  deltaindi     <- array(0,c(p,n,nbiterem))
  tempderiveeas <- matrix(0,p,n)
  dimstatexh    <- 5
  statexh       <- matrix(0,dimstatexh,nbiterem)
  mco           <- matrix(NA,c(n,j))
  mco2          <- array(NA,dim=c(nbiterem,n,j))
  mco3          <- array(NA,dim=c(nbiterem,n,j))
  mcos          <- rep(NA,n)
  mcos2         <- matrix(NA,nbiterem,n)
  STATEXH       <- rep(NA,dimstatexh)
  psisauv       <- array(0,c(n,2,nbiterem))
  
  
  H              <- array(0,c(p,p,n,nbiterem))
  G2             <- array(0,c(p,p,nbiterem))
  tempderiveeas2 <-matrix(0,p,p)
  
  # initial values for the individual parameters
  currentka  <- log(kapop) + rnorm(n,0,sqrt(eta2ka))
  currentcl  <- log(clpop) + rnorm(n,0,sqrt(eta2cl))
  currentpsi <- cbind(exp(currentka),exp(currentcl))
  
  
  ## Start of the EM loop
  for (l in 1:(nbiterem-1)){
    
    ## Simulation step 
    for (k in 1:(n)){
      
      # Parameter ka
      candidatka    <- currentka
      candidatka[k] <- candidatka[k] + rnorm(1,0,sqrt(eta2ka))
      psicandidat   <- cbind(exp(candidatka),exp(currentcl))
      logs          <- -1/2/sigma2*sum((y-model1cptV(psicandidat,id,xidep,vpop))^2)+
        1/2/sigma2*sum((y-model1cptV(currentpsi,id,xidep,vpop))^2)
      logs          <- logs-
        1/2/omega2ka*((candidatka[k]-log(kapop))^2-(currentka[k]-log(kapop))^2)
      u             <- runif(1)
      logu          <- log(u)
      ind           <- (logu<logs)
      currentpsi    <- psicandidat*ind + currentpsi*(1-ind)
      currentka     <- candidatka*ind + currentka*(1-ind)
      
      # Parameter cl
      candidatcl    <- currentcl
      candidatcl[k] <- candidatcl[k] + rnorm(1,0,sqrt(eta2cl))
      psicandidat   <- cbind(exp(currentka),exp(candidatcl))
      logs          <- -1/2/sigma2*sum((y-model1cptV(psicandidat,id,xidep,vpop))^2)+
        1/2/sigma2*sum((y-model1cptV(currentpsi,id,xidep,vpop))^2)
      logs          <- logs -
        1/2/omega2cl*((candidatcl[k]-log(clpop))^2-(currentcl[k]-log(clpop))^2)
      u             <- runif(1)
      logu          <- log(u)
      ind           <- (logu<logs)
      currentpsi    <- psicandidat*ind + currentpsi*(1-ind)
      currentcl     <- candidatcl*ind + currentcl*(1-ind)
      
      # saving simulated data 
      psisauv[,,l]  <- currentpsi
    }
    
    psi <- psisauv[,,l]
    
    ## Parameter estimation update
    
    # estimation of the fixed effect by numerical optimization
    resvpop <- optimize(interval=c(0.001,50),f=floglik,y=y,psi=psisauv[,,1:l],
                        xidep=xidep,id=id,alpha=cumgamma[l,1:l])
    vpop    <- resvpop$minimum
    
    # stochastic approximation of exhaustive statistics and estimation of the 
    # other parameters
    
    mco           <- matrix((y-model1cptV(psi,id,xidep,vpop))^2,n,j,byrow=TRUE)
    mcos          <- apply(mco,1,sum)    
    STATEXH       <- c(apply(log(psi),2,mean), apply(log(psi)^2,2,mean), sum(mcos))
    statexh[,l+1] <- statexh[,l]*(1-gamma[l])+gamma[l]*STATEXH
    
    kapop           <- exp(statexh[1,l+1])
    clpop           <- exp(statexh[2,l+1])
    omega2ka        <- statexh[3,l+1]-statexh[1,l+1]^2
    omega2cl        <- statexh[4,l+1]-statexh[2,l+1]^2
    sigma2          <- statexh[5,l+1]/n/j
    thetaest[,l+1]  <- c(kapop, vpop, clpop, omega2ka, omega2cl, sigma2)
    
    eta2cl <- kRW*omega2cl
    eta2ka <- kRW*omega2ka
    
    
    ## Stochastic approximation of the derivatives of the complete log-likelihood
    
    mco2[l,,]    <- 
      matrix(dVmodel1cpt(psi,id,xidep,vpop)*(y-model1cptV(psi,id,xidep,vpop))/sigma2,
             n,j,byrow=TRUE)
    mcos2[l,]    <- apply(mco2[l,,],1,sum)
    mco3[l,,]    <- 
      matrix((d2Vmodel1cpt(psi,id,xidep,vpop)*(y-model1cptV(psi,id,xidep,vpop))-
                dVmodel1cpt(psi,id,xidep,vpop)^2)/sigma2,n,j,byrow=TRUE)
    
    tempderiveeas[1,] <- (log(psi[,1])-log(kapop))/omega2ka/kapop
    tempderiveeas[3,] <- (log(psi[,2])-log(clpop))/omega2cl/clpop
    tempderiveeas[4,] <- -1/2/omega2ka +1/2/omega2ka^2*(log(psi[,1])-log(kapop))^2
    tempderiveeas[5,] <- -1/2/omega2cl +1/2/omega2cl^2*(log(psi[,2])-log(clpop))^2
    tempderiveeas[6,] <- -j/2/sigma2+apply(mco,1,sum)/2/sigma2^2
    tempderiveeas[2,] <- mcos2[l,]
    
    
    tempderiveeas2[1,1] <-  sum(-log(psi[,1])+log(kapop)-1)/omega2ka/kapop^2
    tempderiveeas2[2,2] <-  sum(mco3[l,,])
    tempderiveeas2[3,3] <-  sum(-log(psi[,2])+log(clpop)-1)/omega2cl/clpop^2
    tempderiveeas2[4,4] <-  n/2/omega2ka^2 -
      1/omega2ka^3*sum((log(psi[,1])-log(kapop))^2)
    tempderiveeas2[5,5] <-  n/2/omega2cl^2 -
      1/omega2cl^3*sum((log(psi[,2])-log(clpop))^2)
    tempderiveeas2[6,6] <-  n*j/2/sigma2^2-
      sum((y-model1cptV(psi,id,xidep,vpop))^2)/sigma2^3
    tempderiveeas2[1,4] <- -sum(log(psi[,1])-log(kapop))/omega2ka^2/kapop
    tempderiveeas2[3,5] <- -sum(log(psi[,2])-log(clpop))/omega2cl^2/clpop
    tempderiveeas2[2,6] <- -sum(mcos2[l,])/(2*sigma2)
    tempderiveeas2[4,1] <- tempderiveeas2[1,4]
    tempderiveeas2[5,3] <- tempderiveeas2[3,5]
    tempderiveeas2[6,2] <- tempderiveeas2[2,6]
    
    
    deltaindi[,,l+1] <- deltaindi[,,l]*(1-gamma[l])+gamma[l]*tempderiveeas
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
  
  res <- list(isco = isco, iobs=iobs, thetaest = thetaest)
  
  
  return(res)
}
