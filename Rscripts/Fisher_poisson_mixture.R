### "Computing an empirical Fisher information matrix estimate in latent variable
### models through stochastic approximation."
### ---------------------------------------------------------------------------

### Numerical study in a Poisson mixture model : statistical properties of 
### two estimators of the Fisher Information matrix 
### (section 4.1 of the article)

## Source useful libraries ----------------------------------------------------
## 
rm(list=ls())

library(ggplot2)
library(cowplot)
library(flexxmix)

## Defining R functions for MCMC estimation of the Fisher Information matrix, 
## computing Iobs and Isco in the Poisson mixture model -----------------------

## Monte-Carlo estimation of the true Fisher information matrix

pm_fisher_mcmc <- function(lambda,alpha,nMC){
  # lambda : vector of Poisson parameter values (as many values as mixture
  # components)
  # alpha  : vector of mixture weights (one less value than components in the
  # mixture, the last weight is deduced from the first)
  # nMC    : 

  z <- rep(0,nMC) # latent variable
  y <- rep(0,nMC) # observed variable
  
  t <- cumsum(c(alpha,1-sum(alpha)))
  
  # data simulation
  for (i in 1:nMC){
    u    <- runif(1)
    z[i] <- 1+sum(u>t)
    y[i] <- rpois(1,lambda[z[i]])
  }
  
  # very large sample based estimation of the FIM
  MCMC_FIM <-  pm_fisher_estimation(y, lambda, alpha)
  MCMC_FIM$Iobs
  MCMC_FIM$Isco
  MCMC_FIM <- (MCMC_FIM$Isco+MCMC_FIM$Iobs)/2
  
  return(MCMC_FIM)
}



## Function for Fisher Information matrix estimation, based on Isco and Iobs, in
##  Poisson mixture models

pm_fisher_estimation <- function(y, lambda, alpha) {
  # y      : vector of observations
  # lambda : vector of mean Poisson parameters for each component of the mixture
  # alpha  : vector of mixture proportions, excluding the proportion of the last
  #  mixture component
  
  K <- length(lambda) # number of components of the mixture
  n <- length(y)      # sample size
  
  deriv1ind  <- matrix(0,2*K-1,n) 
  deriv2     <- matrix(0,2*K-1,2*K-1) 
  covderiv   <- matrix(0,2*K-1,2*K-1)
  
  ## computation of conditional expectation of the first derivatives of the 
  ## complete data log-likelihood
  
  denom <- 0
  for (k in 1:(K-1)){
    denom <- denom + exp(-lambda[k])*lambda[k]^y*alpha[k]
  }
  denom <- denom + exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))
  
  for (k in 1:(K-1)){
    deriv1ind[k,]   <- exp(-lambda[k])*lambda[k]^y*alpha[k]/denom * (y/lambda[k]-1)
    deriv1ind[K+k,] <- exp(-lambda[k])*lambda[k]^y*alpha[k]/denom/alpha[k] + 
      exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(-1/(1-sum(alpha)))
  }
  deriv1ind[K,] <- exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom * (y/lambda[K]-1)
  

  ## computation of conditional expectation of the second derivatives of the 
  ## complete data log-likelihood
  
  
  for (k in 1:(K-1)){
    deriv2[k,k]     <- sum(exp(-lambda[k])*lambda[k]^y*alpha[k]/denom * (-y/lambda[k]^2))
    deriv2[K+k,K+k] <- sum(-exp(-lambda[k])*lambda[k]^y*alpha[k]/denom/alpha[k]^2 
                                  + exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(-1/(1-sum(alpha))^2))
  }
  
  for (k in 1:(K-2)){
    for (l in (k+1):(K-1)){
      deriv2[K+k,K+l] <- sum(exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(-1/(1-sum(alpha))^2))
      deriv2[K+l,K+k] <- deriv2[K+k,K+l] 
    }
  }
  
  deriv2[K,K]<-sum(exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom * (-y/lambda[K]^2))
  
  ## computation of the conditional covariance matrix of the first derivatives of 
  ## the complete data log-likelihood
  
  
  for (k in 1:(K-2)){
    covderiv[k,k] <- sum(exp(-lambda[k])*lambda[k]^y*alpha[k]/denom*(-1+y/lambda[k])^2) 
    covderiv[k+K,k+K] <- sum(exp(-lambda[k])*lambda[k]^y*alpha[k]/denom/alpha[k]^2
                             +exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom/(1-sum(alpha))^2) 
    for (l in (k+1):(K-1)){
      covderiv[k+K,l+K] <- sum(exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom/(1-sum(alpha))^2) 
      covderiv[l+K,k+K] <- covderiv[k+K,l+K]
    } 
    covderiv[k,K+k] <- sum(exp(-lambda[k])*lambda[k]^y*alpha[k]/denom/alpha[k]*(-1+y/lambda[k])) 
    covderiv[k+K,k] <- covderiv[k,K+k]
    
    covderiv[K,K+k] <- sum(exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(-1+y/lambda[K])*(-1)/(1-sum(alpha))) 
    covderiv[K+k,K] <- covderiv[K,K+k]
  }
  
  covderiv[K-1,K-1] <- sum(exp(-lambda[K-1])*lambda[K-1]^y*alpha[K-1]/denom*(-1+y/lambda[K-1])^2) 
  covderiv[2*K-1,2*K-1] <- sum(exp(-lambda[K-1])*lambda[K-1]^y*alpha[K-1]/denom/alpha[K-1]^2+
                                 exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom/(1-sum(alpha))^2) 
  covderiv[K-1,2*K-1] <- sum(exp(-lambda[K-1])*lambda[K-1]^y*alpha[K-1]/denom/alpha[K-1]*(-1+y/lambda[K-1])) 
  covderiv[2*K-1,K-1] <- covderiv[K-1,2*K-1]
  
  covderiv[K,2*K-1] <- sum(exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(-1+y/lambda[K])*(-1)/(1-sum(alpha))) 
  covderiv[2*K-1,K] <- covderiv[K,2*K-1]
  
  covderiv[K,K] <- sum(exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(-1+y/lambda[K])^2) 
  
  ## computation of Isco and Iobs
  Isco <- deriv1ind%*%t(deriv1ind)/n
  Iobs <- deriv1ind%*%t(deriv1ind)/n - deriv2/n - covderiv/n
  
  
  res <- list(Isco = Isco, Iobs = Iobs)
  
  return(res)
}


## Numerical experiment -------------------------------------------------------


nbsim  <- 500        # number of replicates
nMC    <- 100000000  # sample size for Monte-Carlo estimation of the FIM
alpha  <- c(0.3,0.5) # mixture weights of the first K-1 components 
lambda <- c(2,5,9)   # parameter values of the K Poisson distribution of the mixture
n      <- 500        # sample size

## MCMC estimation of the FIM

fisher.mcmc <- pm_fisher_mcmc(lambda,alpha,nMC)

## Computation of Isco and Iobs based on nbsim similated samples

resIobs <- array(NA,dim=c(5,5,nbsim))
resIsco <- array(NA,dim=c(5,5,nbsim))
  
for (j in 1:nbsim){
  
  z <- rep(0,n) # latent variable
  y <- rep(0,n) # observed variable
  
  t <- cumsum(c(alpha,1-sum(alpha)))
  
  for (i in 1:n){
    u    <- runif(1)
    z[i] <- 1+sum(u>t)
    y[i] <- rpois(1,lambda[z[i]])
  }
  
  res <- pm_fisher_estimation(y, lambda, alpha)
  resIobs[,,j] <- res$Iobs
  resIsco[,,j] <- res$Isco
  

}



DataResPoissonMixture <- data.frame(EstF11=c(sqrt(n)*(resIsco[1,1,]-trueFIM[1,1]),sqrt(n)*(resIobs[1,1,]-trueFIM[1,1])),
                         EstF22=c(sqrt(n)*(resIsco[2,2,]-trueFIM[2,2]),sqrt(n)*(resIobs[2,2,]-trueFIM[2,2])),
                         EstF33=c(sqrt(n)*(resIsco[3,3,]-trueFIM[3,3]),sqrt(n)*(resIobs[3,3,]-trueFIM[3,3])),
                         EstF44=c(sqrt(n)*(resIsco[4,4,]-trueFIM[4,4]),sqrt(n)*(resIobs[4,4,]-trueFIM[4,4])),
                         EstF55=c(sqrt(n)*(resIsco[5,5,]-trueFIM[5,5]),sqrt(n)*(resIobs[5,5,]-trueFIM[5,5])),
                         EstF12=c(sqrt(n)*(resIsco[1,2,]-trueFIM[1,2]),sqrt(n)*(resIobs[1,2,]-trueFIM[1,2])),
                         EstF13=c(sqrt(n)*(resIsco[1,3,]-trueFIM[1,3]),sqrt(n)*(resIobs[1,3,]-trueFIM[1,3])),
                         EstF23=c(sqrt(n)*(resIsco[2,3,]-trueFIM[2,3]),sqrt(n)*(resIobs[2,3,]-trueFIM[2,3])),
                         EstF35=c(sqrt(n)*(resIsco[3,5,]-trueFIM[3,5]),sqrt(n)*(resIobs[3,5,]-trueFIM[3,5])),
                         EstF34=c(sqrt(n)*(resIsco[3,4,]-trueFIM[3,4]),sqrt(n)*(resIobs[3,4,]-trueFIM[3,4])),
                         EstF25=c(sqrt(n)*(resIsco[2,5,]-trueFIM[2,5]),sqrt(n)*(resIobs[2,5,]-trueFIM[2,5])),
                         EstF24=c(sqrt(n)*(resIsco[2,4,]-trueFIM[2,4]),sqrt(n)*(resIobs[2,4,]-trueFIM[2,4])),
                         EstF15=c(sqrt(n)*(resIsco[1,5,]-trueFIM[1,5]),sqrt(n)*(resIobs[1,5,]-trueFIM[1,5])),
                         EstF14=c(sqrt(n)*(resIsco[1,4,]-trueFIM[1,4]),sqrt(n)*(resIobs[1,4,]-trueFIM[1,4])),
                         Estimate=c(rep('I N,sco',500),rep('I N,obs',500)))


F11 <- ggplot(DataRes, aes(EstF11, color=Estimate)) +
  geom_density(alpha=.6) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  xlim(-1,1) +
  ggtitle(bquote('('~lambda[1]~','~ lambda[1]~')')) +
  theme(legend.position = c(-0.05, 0.95), plot.title = element_text(size=20,face="bold"))
  #theme(legend.position = "none", plot.title = element_text(size=20,face="bold"))

F22 <- ggplot(DataRes, aes(EstF22, color=Estimate)) +
  geom_density(alpha=.6) +
  scale_fill_manual(values = c("#984EA3",'#E69F00'))  +
  xlab("") +
  ylab("") +
  xlim(-0.5,0.5) +
  ggtitle(bquote('('~lambda[2]~','~ lambda[2]~')')) +
  theme(legend.position = "none", plot.title = element_text(size=20,face="bold"))

F25 <- ggplot(DataRes, aes(EstF25, color=Estimate)) +
  geom_density(alpha=.6) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  xlim(-1.5,1.5) +
  ggtitle(bquote('('~lambda[2]~','~ alpha[2]~')')) +
  theme(legend.position = "none", plot.title = element_text(size=20,face="bold"))


plot_grid(F11, F22, F25, ncol = 3, nrow = 1)
