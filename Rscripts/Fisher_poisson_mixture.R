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
library(flexmix)

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
    deriv1ind[k,]   <- exp(-lambda[k])*lambda[k]^y*alpha[k]/denom*(y/lambda[k]-1)
    deriv1ind[K+k,] <- exp(-lambda[k])*lambda[k]^y/denom - 
      exp(-lambda[K])*lambda[K]^y/denom
  }
  deriv1ind[K,] <- exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(y/lambda[K]-1)
  
  
  ## computation of conditional expectation of the second derivatives of the 
  ## complete data log-likelihood
  
  
  for (k in 1:(K-1)){
    deriv2[k,k]     <- sum(exp(-lambda[k])*lambda[k]^y*alpha[k]/denom*(-y/lambda[k]^2))
    deriv2[K+k,K+k] <- sum(-exp(-lambda[k])*lambda[k]^y/denom/alpha[k] -
                           exp(-lambda[K])*lambda[K]^y/denom*(1/(1-sum(alpha))))
  }
  
  for (k in 1:(K-2)){
    for (l in (k+1):(K-1)){
      deriv2[K+k,K+l] <- - sum(exp(-lambda[K])*lambda[K]^y/denom*(1/(1-sum(alpha))))
      deriv2[K+l,K+k] <- deriv2[K+k,K+l] 
    }
  }
  
  deriv2[K,K]<-sum(exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(-y/lambda[K]^2))
  
  ## computation of the conditional covariance matrix of the first derivatives of 
  ## the complete data log-likelihood
  
  
  for (k in 1:(K-2)){
    covderiv[k,k] <- sum(exp(-lambda[k])*lambda[k]^y*alpha[k]/denom*(-1+y/lambda[k])^2) 
    covderiv[k+K,k+K] <- sum(exp(-lambda[k])*lambda[k]^y/denom/alpha[k] +
                             exp(-lambda[K])*lambda[K]^y/denom/(1-sum(alpha))) 
    for (l in (k+1):(K-1)){
      covderiv[k+K,l+K] <- sum(exp(-lambda[K])*lambda[K]^y/denom/(1-sum(alpha))) 
      covderiv[l+K,k+K] <- covderiv[k+K,l+K]
    } 
    covderiv[k,K+k] <- sum(exp(-lambda[k])*lambda[k]^y/denom*(-1+y/lambda[k])) 
    covderiv[k+K,k] <- covderiv[k,K+k]
    
    covderiv[K,K+k] <- sum(exp(-lambda[K])*lambda[K]^y/denom*(-1+y/lambda[K])*(-1)) 
    covderiv[K+k,K] <- covderiv[K,K+k]
  }
  
  covderiv[K-1,K-1] <- sum(exp(-lambda[K-1])*lambda[K-1]^y*alpha[K-1]/denom*(-1+y/lambda[K-1])^2) 
  covderiv[2*K-1,2*K-1] <- sum(exp(-lambda[K-1])*lambda[K-1]^y/denom/alpha[K-1]+
                                 exp(-lambda[K])*lambda[K]^y/denom/(1-sum(alpha))) 
  covderiv[K-1,2*K-1] <- sum(exp(-lambda[K-1])*lambda[K-1]^y/denom*(-1+y/lambda[K-1])) 
  covderiv[2*K-1,K-1] <- covderiv[K-1,2*K-1]
  
  covderiv[K,2*K-1] <- sum(exp(-lambda[K])*lambda[K]^y/denom*(-1+y/lambda[K])*(-1)/(1-sum(alpha))) 
  covderiv[2*K-1,K] <- covderiv[K,2*K-1]
  
  covderiv[K,K] <- sum(exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(-1+y/lambda[K])^2) 
  
  ## computation of Isco and Iobs
  Isco <- deriv1ind%*%t(deriv1ind)/n
  Iobs <- deriv1ind%*%t(deriv1ind)/n - deriv2/n - covderiv/n
  
  
  res <- list(Isco = Isco, Iobs = Iobs)
  
  return(res)
}


## Numerical experiment -------------------------------------------------------


nbsim  <- 100 #500        # number of replicates
nMC    <- 100000000  # sample size for Monte-Carlo estimation of the FIM
alpha  <- c(0.3,0.5) # mixture weights of the first K-1 components 
lambda <- c(2,5,9)   # parameter values of the K Poisson distribution of the mixture
n      <- 500        # sample size

theta.true <- matrix(c(lambda,alpha),ncol=1)

## MCMC estimation of the FIM

fisher.mcmc <- pm_fisher_mcmc(lambda,alpha,nMC)
inv.sqrt.fisher.mcmc <- solve(chol(fisher.mcmc))
## Cholesky de fisher.mcmc à caluler une seule fois, en dehors de la boucle

## Computation of Isco and Iobs based on nbsim simulated samples

Iobs.true.theta <- array(NA,dim=c(5,5,nbsim))
Isco.true.theta <- array(NA,dim=c(5,5,nbsim))
Iobs.theta.est <- array(NA,dim=c(5,5,nbsim))
Isco.theta.est <- array(NA,dim=c(5,5,nbsim))  

est.lambda <- matrix(NA,3,nbsim)
est.alpha <- matrix(NA,2,nbsim)

## counters for the estimation of the coverage rates
coverage.iobs.theta.est <- 0
coverage.isco.theta.est <- 0

coverage.iobs.theta.true <- 0
coverage.isco.theta.true <- 0

coverage.true.fisher <- 0

failure.cov.Isco.theta.est <- 0

## TBA

failed <- 0
det.Iobs.theta.true <- 0
det.Iobs.theta.est <- 0
det.Isco.theta.true <- 0
det.Isco.theta.est <- 0

j <- 1

while (j <= nbsim){
  
  ## Data simulation
  
  z <- rep(0,n) # latent variable
  y <- rep(0,n) # observed variable
  
  t <- cumsum(c(alpha,1-sum(alpha)))
  
  for (i in 1:n){
    u    <- runif(1)
    z[i] <- 1+sum(u>t)
    y[i] <- rpois(1,lambda[z[i]])
  }
  
  ## Computation of Isco and Iobs in the true parameter value
  
  res.true.theta <- pm_fisher_estimation(y, lambda, alpha)
  Iobs.true.theta[,,j] <- res.true.theta$Iobs
  Isco.true.theta[,,j] <- res.true.theta$Isco
  
  ## Computation of Isco and Iobs in the MLE value of the parameter
  
  ### Maximum likelihood estimation
  
  em.est <- flexmix(y ~ 1, k = 3, model = FLXMCmvpois(),
                    control = list(tolerance = 1e-15, iter.max = 1000,
                                   nrep=10))
  if ((length(parameters(em.est))==3) & (sum(diff(parameters(em.est))>0.1)==2)){
    ## Both conditions are to avoid identifiability issues
    est.lambda[,j] <- parameters(em.est)
    ord <- order(est.lambda[,j]) ## Fix label switching issues
    est.lambda[,j] <- est.lambda[ord,j]
    est.alpha[,j] <- prior(em.est)[ord[-3]]
    
    ### Estimation of the FIM
    
    res.theta.est <- pm_fisher_estimation(y, est.lambda[,j], est.alpha[,j])
    Iobs.theta.est[,,j] <- res.theta.est$Iobs
    Isco.theta.est[,,j] <- res.theta.est$Isco
    
    
    ## Coverage rates of confidence ellipsoids ...
    theta.est <- matrix(c(est.lambda[,j],est.alpha[,j]),ncol=1)
    unit.vect <- matrix(rep(1,length(theta.est)),ncol=1)
    
    ## ... based on the true Fisher information matrix
    
    conf.inf.fisher.true <- theta.est -
      1.96/sqrt(n)*inv.sqrt.fisher.mcmc%*%unit.vect
    conf.sup.fisher.true <- theta.est +
      1.96/sqrt(n)*inv.sqrt.fisher.mcmc%*%unit.vect
    
    coverage.true.fisher <- coverage.true.fisher +
      (theta.true<=conf.sup.fisher.true)*(theta.true>=conf.inf.fisher.true)
    
    
    ## ... based on Isco computed in the MLE value of the parameters
    
    #tryCatch(
    #  {
    if (det(Isco.theta.est[,,j])>0){
      chol.Isco.t.e <- chol(Isco.theta.est[,,j],p=1)
      conf.inf.isco.theta.est <- theta.est -
        1.96/sqrt(n)*solve(chol.Isco.t.e)%*%unit.vect
      conf.sup.isco.theta.est <- theta.est +
        1.96/sqrt(n)*solve(chol.Isco.t.e)%*%unit.vect
      
      coverage.isco.theta.est <- coverage.isco.theta.est +
        (theta.true<=conf.sup.isco.theta.est)*(theta.true>=conf.inf.isco.theta.est)
    } else{
      det.Isco.theta.est <- det.Isco.theta.est+1
    }
    #},
    #finally = {
    #  failure.cov.Isco.theta.est <- failure.cov.Isco.theta.est +1
    #}
    #)
    
    ## ... based on Isco computed in the true parameter values
    
    if (det(Isco.true.theta[,,j])>0){
      chol.Isco.t.t <- chol(Isco.true.theta[,,j],p=1)
      conf.inf.isco.theta.true <- theta.est -
        1.96/sqrt(n)*solve(chol.Isco.t.t)%*%unit.vect
      conf.sup.isco.theta.true <- theta.est +
        1.96/sqrt(n)*solve(chol.Isco.t.t)%*%unit.vect
      
      coverage.isco.theta.true <- coverage.isco.theta.true +
        (theta.true<=conf.sup.isco.theta.true)*(theta.true>=conf.inf.isco.theta.true)
    } else{
      det.Isco.theta.true <- det.Isco.theta.true+1
    }
    ## ... based on Iobs computed in the MLE value of the parameters
    
    if (det(Iobs.theta.est[,,j])>0){
      chol.Iobs.t.e <- chol(Iobs.theta.est[,,j],p=1)
      conf.inf.iobs.theta.est <- theta.est -
        1.96/sqrt(n)*solve(chol.Iobs.t.e)%*%unit.vect
      conf.sup.iobs.theta.est <- theta.est +
        1.96/sqrt(n)*solve(chol.Iobs.t.e)%*%unit.vect
      
      coverage.iobs.theta.est <- coverage.iobs.theta.est +
        (theta.true<=conf.sup.iobs.theta.est)*(theta.true>=conf.inf.iobs.theta.est)
    } else{
      det.Iobs.theta.est <- det.Iobs.theta.est+1
    }
    ## ... based on Iobs computed in the true parameter values
    
    if (det(Iobs.true.theta[,,j])>0){
      chol.Iobs.t.t <- chol(Iobs.true.theta[,,j],p=1)
      conf.inf.iobs.theta.true <- theta.est -
        1.96/sqrt(n)*solve(chol.Iobs.t.t)%*%unit.vect
      conf.sup.iobs.theta.true <- theta.est +
        1.96/sqrt(n)*solve(chol.Iobs.t.t)%*%unit.vect
      
      coverage.iobs.theta.true <- coverage.iobs.theta.true +
        (theta.true<=conf.sup.iobs.theta.true)*(theta.true>=conf.inf.iobs.theta.true)
    } else{
      det.Iobs.theta.true <- det.Iobs.theta.true+1
    }
    
    j <- j+1
  } else{
    failed <- failed +1
  }
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
