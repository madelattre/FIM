## R script for studying the relevance of the SAEM algorithm out of the curved 
## exponential family 
## -----------------------------------------------------------------------------


## 1- Data simulation

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

## Estimation

nbiterem     <- 3000 # total number of iterations
nbiterburnin <- 1000 # number of burnin iterations

nbsim <- 500  # number of simulated datasets  

thetaest <- matrix(0,6,nbsim)
isco.est <- array(0,dim=c(6,6,nbsim))
iobs.est <- array(0,dim=c(6,6,nbsim))

for (kk in 1:nbsim){
  
  set.seed(kk*2500+10)
  
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
  
  datasim <- data.frame(y=y,dose=rep(dose,n*j),time=rep(times,n),
                        subject=kronecker(1:n, rep(1,j)))
  
  ## Estimation
  theta0 <- list(vpop=vpop*runif(1,0.95,1.05),kapop=kapop*runif(1,0.8,1.2),
                 clpop=clpop*runif(1,0.8,1.2),omega2ka=omega2ka*runif(1,0.4,2),
                 omega2cl=omega2cl*runif(1,0.4,2),sigma2=sigma2*runif(1,0.4,2))
  
  res <- saem_non_exp(datasim, nbiterem, nbiterburnin, theta0)
  
  thetaest[,kk]  <- res$thetaest[,nbiterem]
  isco.est[,,kk] <- res$isco[,,nbiterem]
  iobs.est[,,kk] <- res$iobs[,,nbiterem]
  
  filename <- paste("Rfiles/ResNLMEnonexponential.Rdata",sep="")
  resNLMEnonExp <- list(thetaest=thetaest,isco=isco.est,iobs=iobs.est)
  save(resNLMEnonExp,file=filename)
}
