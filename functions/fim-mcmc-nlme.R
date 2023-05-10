## R function for computing Monte-carlo estimation of the FIM in the PK nonlinear
## mixed-effects model with three random effects belonging to the curved 
## exponential family

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