## R script for studying the properties of the SAEM algorithm in the curved 
## exponential family when the number of iterations grow. 
## -----------------------------------------------------------------------------

## 1- Data simulation

# Sample characteristics

set.seed(3005)

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

datasim <- data.frame(y=y,dose=rep(dose,n*j),time=rep(times,n),
                      subject=kronecker(1:n, rep(1,j)))


## 2- Numerical experiment

## a- Evaluation of both estimators of the FIM using the saem algorithm

nbsim <- 500

# Algorithmic settings
nbiterem <- 3000
nbiterburnin <- 1000

# Saving the nbsim results
iscoarray <- array(0,dim=c(nbsim,7,7,nbiterem))
iobsarray <- array(0,dim=c(nbsim,7,7,nbiterem))
thetaest  <- matrix(NA,7,nbsim)

for (k in 1:nbsim){
  
  set.seed(k*100+10)
  
  theta0         <- list(vpop=vpop*runif(1,0.8,1.2),
                         kapop=kapop*runif(1,0.8,1.2),
                         clpop=clpop*runif(1,0.8,1.2),
                         omega2v=omega2v*runif(1,0.4,2),
                         omega2ka=omega2ka*runif(1,0.4,2),
                         omega2cl=omega2cl*runif(1,0.4,2),
                         sigma2=sigma2*runif(1,0.4,2))
  res            <- saem(datasim, nbiterem, nbiterburnin, theta0)
  iscoarray[k,,,]<- res$isco
  iobsarray[k,,,]<- res$iobs
  thetaest[,k]   <- res$thetaest[,nbiterem]
  
}

# b- Monte-Carlo evaluation of both estimates
# These Monte-Carlo estimations are considered as the targets for the estimates
# computed using the stochastic approximation algorithm

nbMC <- 10000
nbMCburnin <- 5000

tm <- rowMeans(thetaest)

thetaMean <- list(kapop=tm[1],vpop=tm[2],clpop=tm[3],omega2ka=tm[4],
                  omega2v=tm[5],omega2cl=tm[6],sigma2=tm[7])

FisherMC <- FIM_mc(datasim, nbMC, nbMCburnin, thetaMean)

iscoMC <- FisherMC$iscoMC
iobsMC <- FisherMC$iobsMC

# Evaluation of the mean relative bias and of the mean relative standard errors
# per iteration.

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

MbiasIsco <- apply(biasIsco[,,,(nbiterburnin+1):nbiterem],c(1,2,4),mean)
MbiasIobs <- apply(biasIobs[,,,(nbiterburnin+1):nbiterem],c(1,2,4),mean)

MsdIsco <- apply(rsdIsco[,,,(nbiterburnin+1):nbiterem],c(1,2,4),mean)
MsdIobs <- apply(rsdIobs[,,,(nbiterburnin+1):nbiterem],c(1,2,4),mean)


save(MbiasIsco,file='Rfiles/ResNLMEexponentialBiasIsco.Rdata')
save(MbiasIobs,file='Rfiles/ResNLMEexponentialBiasIobs.Rdata')
save(MsdIsco,file='Rfiles/ResNLMEexponentialSdIsco.Rdata')
save(MsdIobs,file='Rfiles/ResNLMEexponentialSdIobs.Rdata')

load('Rfiles/ResNLMEexponentialBiasIsco.Rdata')
load('Rfiles/ResNLMEexponentialBiasIobs.Rdata')
load('Rfiles/ResNLMEexponentialSdIsco.Rdata')
load('Rfiles/ResNLMEexponentialSdIobs.Rdata')

MbiasIobs <- MbiasIobs[,,seq(1,2000,10)]
MbiasIsco <- MbiasIsco[,,seq(1,2000,10)]
MsdIobs <- MsdIobs[,,seq(1,2000,10)]
MsdIsco <- MsdIsco[,,seq(1,2000,10)]

save(MbiasIsco,file='Rfiles/ResNLMEexponentialBiasIsco.Rdata')
save(MbiasIobs,file='Rfiles/ResNLMEexponentialBiasIobs.Rdata')
save(MsdIsco,file='Rfiles/ResNLMEexponentialSdIsco.Rdata')
save(MsdIobs,file='Rfiles/ResNLMEexponentialSdIobs.Rdata')
