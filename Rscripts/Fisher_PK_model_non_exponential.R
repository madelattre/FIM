## Nonlinear PK mixed-effects model not belonging to the curved exponential family 
## Properties of the stochastic approximation algorithm for the evaluation of the two estimators of the Fisher information matrix
## (Section 4.2.1)

rm(list=ls())

## 1-Estimation functions

# SAEM algorithm

saem <- function(data, nbiterem, nbiterburnin, theta0, kRW=0.5) {
  
  # data         : dataset 
  # nbiterem     : total number of iterations of the SAEM algorithm
  # nbiterburnin : number of burn-in iterations of the algorithm
  # theta0       : initial parameter values
  # kRW          : coefficient used to adjust the variance of the proposal kernel of the MCMC procedure
  
  # Important functions used within each iteration of the SAEM algorithm
  
  ## conditional mean of the observations given the random parameters psi, fixed effect V,
  ## and structural covariates, in this case dose and time
  
  model1cpt<-function(psi,id,xidep,V) {
    dose  <- xidep[,1]
    tim   <- xidep[,2]
    ka    <- psi[id,1]
    CL    <- psi[id,2]
    k     <- CL/V
    ypred <- dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
    return(ypred)
  }
  
  # Derivative of the model function with respect to V
  
  dVmodel1cpt<-function(psi,id,xidep,V) {
    dose  <- xidep[,1]
    tim   <- xidep[,2]
    ka    <- psi[id,1]
    CL    <- psi[id,2]
    ypred <- -dose*ka^2/(V*ka-CL)^2*(exp(-CL/V*tim)-exp(-ka*tim))+dose*ka/(V*ka-CL)*CL/V^2*tim*exp(-CL/V*tim)
    return(ypred)
  }
  
  # Q quantity 
  
  floglik <- function(v,y,psi,xidep,id,alpha){
    l <- length(alpha)
    psi <- array(psi,dim=c(n,2,l))
    value <- 0
    for (ll in 1:l){
      moyij <- model1cpt(psi[,,ll],id,xidep,v)
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
  mcos          <- rep(NA,n)
  mcos2         <- matrix(NA,nbiterem,n)
  STATEXH       <- rep(NA,dimstatexh)
  psisauv       <- array(0,c(n,2,nbiterem))

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
      logs          <- -1/2/sigma2*sum((y-model1cpt(psicandidat,id,xidep,vpop))^2)+1/2/sigma2*sum((y-model1cpt(currentpsi,id,xidep,vpop))^2)
      logs          <- logs-1/2/omega2ka*((candidatka[k]-log(kapop))^2-(currentka[k]-log(kapop))^2)
      u             <- runif(1)
      logu          <- log(u)
      ind           <- (logu<logs)
      currentpsi    <- psicandidat*ind + currentpsi*(1-ind)
      currentka     <- candidatka*ind + currentka*(1-ind)
    
      # Parameter cl
      candidatcl    <- currentcl
      candidatcl[k] <- candidatcl[k] + rnorm(1,0,sqrt(eta2cl))
      psicandidat   <- cbind(exp(currentka),exp(candidatcl))
      logs          <- -1/2/sigma2*sum((y-model1cpt(psicandidat,id,xidep,vpop))^2)+1/2/sigma2*sum((y-model1cpt(currentpsi,id,xidep,vpop))^2)
      logs          <- logs -1/2/omega2cl*((candidatcl[k]-log(clpop))^2-(currentcl[k]-log(clpop))^2)
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
    resvpop <- optimize(interval=c(0.001,50),f=floglik,y=y,psi=psisauv[,,1:l],xidep=xidep,id=id,alpha=cumgamma[l,1:l])
    vpop    <- resvpop$minimum
    
    # stochastic approximation of exhaustive statistics and estimation of the other parameters
    
    mco           <- matrix((y-model1cpt(psi,id,xidep,vpop))^2,n,j,byrow=TRUE)
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
    
    mco2[l,,]    <- matrix(dVmodel1cpt(psi,id,xidep,vpop)*(y-model1cpt(psi,id,xidep,vpop))/sigma2,n,j,byrow=TRUE)
    mcos2[l,]    <- apply(mco2[l,,],1,sum)
    
    tempderiveeas[1,] <- (log(psi[,1])-log(kapop))/omega2ka/kapop
    tempderiveeas[3,] <- (log(psi[,2])-log(clpop))/omega2cl/clpop
    tempderiveeas[4,] <- -1/2/omega2ka +1/2/omega2ka^2*(log(psi[,1])-log(kapop))^2
    tempderiveeas[5,] <- -1/2/omega2cl +1/2/omega2cl^2*(log(psi[,2])-log(clpop))^2
    tempderiveeas[6,] <- -j/2/sigma2+apply(mco,1,sum)/2/sigma2^2
    tempderiveeas[2,] <- mcos2[l,]
    
    deltaindi[,,l+1] <- deltaindi[,,l]*(1-gamma[l])+gamma[l]*tempderiveeas 
  }
  ## End of the em loop
  
  ## Computation of the FIM estimators based on the score
  
  isco <- array(0,c(p,p,nbiterem))
  
  for (l in 1:nbiterem){
    isco[,,l] <- deltaindi[,,l]%*%t(deltaindi[,,l])/n
  }
  
  
  res <- list(isco = isco, thetaest = thetaest)
  
  
  return(res)
}

##########################
##########################

## 2-Data simulation

## Model function 
## returns the conditional mean of the observations given the random parameters psi, time and dose

model1cptsim<-function(psi,id,tim,dose) {
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  k<-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

## Sample characteristics 

n     <- 100                             # number of subjects
times <- c(0.25,0.5,1,2,3.5,5,7,9,12,24) # observation times
j     <- length(times)                   # number of observations per subject
dose  <- 320                             # dose

## True parameter values
vpop     <- 31
kapop    <- 1.6
clpop    <- 2.8
omega2ka <- 0.40
omega2cl <- 0.40 
sigma2   <- 0.75

## 2.1 Estimation and convergence graphs on one simulated dataset

## Simulation of individual parameters
vi  <- rep(vpop,n)
kai <- exp(rnorm(n,log(kapop),sd=sqrt(omega2ka)))
cli <- exp(rnorm(n,log(clpop),sd=sqrt(omega2cl)))

## Simulation of the observations
ypred <- c()
for (k in 1:n){
  ypred <- c(ypred,model1cptsim(cbind(kai,vi,cli),k,times,dose))
}

y <- ypred + rnorm(n*j,0,sd=sqrt(sigma2))

datasim <- data.frame(y=y,dose=rep(dose,n*j),time=rep(times,n),subject=kronecker(1:n, rep(1,j)))

## Estimation

# Algorithmic settings
nbiterem     <- 3000
nbiterburnin <- 1500
theta0       <- list(vpop=(vpop-5)*runif(1,0.8,1.2),kapop=kapop*runif(1,0.8,1.2),clpop=clpop*runif(1,0.8,1.2),
                 omega2ka=omega2ka*runif(1,0.4,2),omega2cl=omega2cl*runif(1,0.4,2), sigma2=sigma2*runif(1,0.4,2))

# Parameter and FIM estimation 
res <- saem(datasim, nbiterem, nbiterburnin, theta0)


# Convergence graphs

library(ggplot2)
library(cowplot)
library(devtools)
library(RColorBrewer)

DataResEst <- data.frame(
  ka      = res$thetaest[1,],
  V        = res$thetaest[2,],  
  Cl       = res$thetaest[3,],
  o2ka     = res$thetaest[4,],
  Iscoka   = res$isco[1,1,],
  Iscov    = res$isco[2,2,],
  Iscocl   = res$isco[3,3,],
  Iscoo2ka = res$isco[4,4,],
  Iter=seq(1,nbiterem))

kaconv <- ggplot(DataResEst, aes(y=ka, x=Iter)) +
  geom_line(size=0.5) +
  xlab("") +
  ylab("") +
  ggtitle(bquote(ka)) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=45))

vconv <- ggplot(DataResEst, aes(y=V, x=Iter)) +
  geom_line(size=0.5) +
  xlab("") +
  ylab("") +
  ggtitle(bquote(V)) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=45))

clconv <- ggplot(DataResEst, aes(y=Cl, x=Iter)) +
  geom_line(size=0.5) +
  xlab("") +
  ylab("") +
  ylim(2.5,3) +
  ggtitle(bquote(Cl)) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=45))

o2kaconv <- ggplot(DataResEst, aes(y=o2ka, x=Iter)) +
  geom_line(size=0.5) +
  xlab("") +
  ylab("") +
  ggtitle(bquote(omega[ka]^2)) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=45))

iscokaconv <- ggplot(DataResEst, aes(y=Iscoka, x=Iter)) +
  geom_line(size=0.5) +
  xlab("") +
  ylab("") +
  ggtitle(bquote(I[n-sco](ka,ka))) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=45))

iscovconv <- ggplot(DataResEst, aes(y=Iscov, x=Iter)) +
  geom_line(size=0.5) +
  xlab("") +
  ylab("") +
  ggtitle(bquote(I[n-sco](V,V))) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=45))

iscoclconv <- ggplot(DataResEst, aes(y=Iscocl, x=Iter)) +
  geom_line(size=0.5) +
  xlab("") +
  ylab("") +
  ggtitle(bquote(I[n-sco](Cl,Cl))) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=45))

iscoo2kaconv <- ggplot(DataResEst, aes(y=Iscoo2ka, x=Iter)) +
  geom_line(size=0.5) +
  xlab("") +
  ylab("") +
  ggtitle(bquote(I[n-sco](omega[ka]^2,omega[ka]^2))) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=45))

plot_grid(kaconv, vconv, clconv, o2kaconv, iscokaconv, iscovconv, iscoclconv, iscoo2kaconv, ncol = 4, nrow = 2)


## 2.2 Empirical covering rates from nbsim simulated datasets

nbiterem     <- 3000 # total number of iterations
nbiterburnin <- 1000 # number of burnin iterations

nbsim <- 500  # number of simulated datasets  
rate  <- 0.95 # nominal rate

covering.isco <- rep(0,5)

for (kk in 1:nbsim){
  print(kk)
  
  ## Simulation of individual parameters
  vi  <- rep(vpop,n)
  kai <- exp(rnorm(n,log(kapop),sd=sqrt(omega2ka)))
  cli <- exp(rnorm(n,log(clpop),sd=sqrt(omega2cl)))
  
  ## Simulation of the observations
  ypred <- c()
  for (k in 1:n){
    ypred <- c(ypred,model1cptsim(cbind(kai,vi,cli),k,times,dose))
  }
  
  y <- ypred + rnorm(n*j,0,sd=sqrt(sigma2))
  
  datasim <- data.frame(y=y,dose=rep(dose,n*j),time=rep(times,n),subject=kronecker(1:n, rep(1,j)))
  
  ## Estimation
  
  res <- saem(datasim, nbiterem, nbiterburnin, theta0)
  
  thetaest <- res$thetaest[,nbiterem]
  isco <- res$isco[,,nbiterem]
  
  ## Confidence intervals 
  
  ICinfisco <- thetaest - qnorm(1-(1-rate)/2,0,1)*sqrt(diag(solve(isco))/n)
  ICsupisco <- thetaest + qnorm(1-(1-rate)/2,0,1)*sqrt(diag(solve(isco))/n)
  
  covering.isco <- covering.isco+as.numeric((theta.true<=ICsupisco)&(theta.true>=ICinfisco))

}

covering.isco/nbsim
