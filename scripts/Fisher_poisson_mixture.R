### "Computing an empirical Fisher information matrix estimate in latent variable
### models through stochastic approximation."
### ---------------------------------------------------------------------------

### Numerical study in a Poisson mixture model : statistical properties of 
### two estimators of the Fisher Information matrix 
### (section 4.1 of the article)

## Source useful libraries ----------------------------------------------------

rm(list=ls())

library(ggplot2)
library(cowplot)


## Defining useful R functions ------------------------------------------------




## Simulation of observations according to a mixture of Poisson distributions

data.sim <- function(n,lambda,alpha){
  # n      : number of observations to be simulated
  # lambda : vector of K Poisson parameters in ascending order
  # alpha  : vector of mixture proportions (first K-1 components)
  
  z <- rep(0,n) # latent variable
  y <- rep(0,n) # observed variable
  
  t <- cumsum(c(alpha,1-sum(alpha)))
  
  for (i in 1:n){
    u    <- runif(1)
    z[i] <- 1+sum(u>t)
    y[i] <- rpois(1,lambda[z[i]])
  }
  
  return(y)
}

## Function for Fisher Information matrix estimation, 
## computing Isco and Iobs

pm_fisher_estimation <- function(y, lambda, alpha) {
  # y      : vector of observations
  # lambda : vector of K Poisson parameters for each component of the mixture
  #         (in ascending order)
  # alpha  : vector of (K-1) mixture proportions
  
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
  
  covderiv[K,2*K-1] <- sum(exp(-lambda[K])*lambda[K]^y/denom*(-1+y/lambda[K])*(-1)) 
  covderiv[2*K-1,K] <- covderiv[K,2*K-1]
  
  covderiv[K,K] <- sum(exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(-1+y/lambda[K])^2) 
  
  ## computation of Isco and Iobs
  
  Isco <- deriv1ind%*%t(deriv1ind)/n
  Iobs <- deriv1ind%*%t(deriv1ind)/n - deriv2/n - covderiv/n
  
  
  res <- list(Isco = Isco, Iobs = Iobs)
  
  return(res)
}


## Monte-Carlo estimation of the true Fisher information matrix

pm_fisher_mcmc <- function(lambda,alpha,nMC){
  # lambda : vector of K Poisson parameter values (in ascending order)
  # alpha  : vector of K-1 mixture weights (one less value than components in the
  # mixture, the last weight is deduced from the first)
  # nMC    : size of the Monte-Carlo sample
  
  z <- rep(0,nMC) # latent variable
  y <- rep(0,nMC) # observed variable
  
  t <- cumsum(c(alpha,1-sum(alpha)))
  
  # data simulation
  y <- data.sim(nMC,lambda,alpha)
  
  # very large sample based estimation of the FIM
  MCMC_FIM <-  pm_fisher_estimation(y, lambda, alpha)
  MCMC_FIM$Iobs
  MCMC_FIM$Isco
  MCMC_FIM <- (MCMC_FIM$Isco+MCMC_FIM$Iobs)/2
  
  return(MCMC_FIM)
}


## Intermediate function for calculating empirical coverage rates and checking 
## the positive definiteness of Fisher estimates

coverage <- function(theta.true,theta.est,FIM,rate,n){
  # theta.true: vector of true parameter values in the following order : K 
  # Poisson parameters in ascending order, K-1 first mixture proportions
  # theta.est: vector of estimated parameter values in the following order : K 
  # Poisson parameters in ascending order, K-1 first mixture proportions
  # FIM: Fisher Information matrix or its estimate
  # rate: coverage rate
  # n: sample size
  
  DP <- 1
  
  quantile <- qnorm((1+rate)/2)
  
  conf.inf <- theta.est - quantile/sqrt(n)*sqrt(diag(solve(FIM)))
  conf.sup <- theta.est + quantile/sqrt(n)*sqrt(diag(solve(FIM)))
  
  coverage <- (theta.true<=conf.sup)*(theta.true>=conf.inf)
  
  if (sum((eigen(FIM)$values)<=0)>0){
    DP <- 0
  }
  
  res <- list(coverage=coverage,DP=DP)
  
  return(res)
  
}



## Numerical experiment -------------------------------------------------------


nbsim  <- 1000       # number of replicates
nMC    <- 1000000    # sample size for Monte-Carlo estimation of the FIM
alpha  <- c(0.3,0.5) # mixture weights of the first K-1 components 
lambda <- c(2,5,9)   # parameter values of the K Poisson distribution of the mixture
n      <- 500        # sample size

rate <- 0.95

theta.true <- matrix(c(lambda,alpha),ncol=1)

## MCMC estimation of the FIM
fisher.mcmc <- pm_fisher_mcmc(lambda,alpha,nMC)

## Objects for storing results
### Estimation of the FIM while knowing the true parameter values or while
### estimating the parameter values
Iobs.true.theta <- array(NA,dim=c(5,5,nbsim))
Isco.true.theta <- array(NA,dim=c(5,5,nbsim))
Iobs.theta.est  <- array(NA,dim=c(5,5,nbsim))
Isco.theta.est  <- array(NA,dim=c(5,5,nbsim))  
### Parameter estimates
est.lambda <- matrix(NA,3,nbsim)
est.alpha  <- matrix(NA,2,nbsim)
### To derive the empirical coverage rate
coverage.fisher.mcmc <- matrix(0,nbsim,5)
coverage.iobstrue    <- matrix(0,nbsim,5)
coverage.iobsest     <- matrix(0,nbsim,5)
coverage.iscotrue    <- matrix(0,nbsim,5)
coverage.iscoest     <- matrix(0,nbsim,5)
### Positive definiteness of the FIM estimates
DP.iobstrue <- rep(0,nbsim)
DP.iobsest  <- rep(0,nbsim)
DP.iscotrue <- rep(0,nbsim)
DP.iscoest  <- rep(0,nbsim)



for (j in 1:nbsim){
  
  ## Data simulation
  
  y <- data.sim(n,lambda,alpha)
  
  ## Computation of Isco and Iobs in the true parameter value
  
  res.true.theta <- pm_fisher_estimation(y, lambda, alpha)
  Iobs.true.theta[,,j] <- res.true.theta$Iobs
  Isco.true.theta[,,j] <- res.true.theta$Isco
  
  ## Computation of Isco and Iobs in the MLE value of the parameter
  
  em.est <- em_poisson_mixture(y,3)
  est.lambda[,j] <- em.est[[1]]
  est.alpha[,j] <- em.est[[2]]#[ord[-3]]
  
  res.theta.est <- pm_fisher_estimation(y, est.lambda[,j], est.alpha[,j])
  Iobs.theta.est[,,j] <- res.theta.est$Iobs
  Isco.theta.est[,,j] <- res.theta.est$Isco
  
  
  if (n>=500){
    ## Coverage rates and positive definiteness of the FIM estimates
    
    cov0 <- coverage(theta.true,c(est.lambda[,j],est.alpha[,j]),fisher.mcmc,rate,n)
    coverage.fisher.mcmc[j,] <- cov0$coverage
    
    cov1 <- coverage(theta.true,c(est.lambda[,j],est.alpha[,j]),
                     Isco.true.theta[,,j],rate,n)
    coverage.iscotrue[j,] <- cov1$coverage
    DP.iscotrue[j] <- cov1$DP
    
    cov2 <- coverage(theta.true,c(est.lambda[,j],est.alpha[,j]),
                     Isco.theta.est[,,j],rate,n)
    coverage.iscoest[j,] <- cov2$coverage
    DP.iscoest[j] <- cov2$DP
    
    cov3 <- coverage(theta.true,c(est.lambda[,j],est.alpha[,j]),
                     Iobs.true.theta[,,j],rate,n)
    coverage.iobstrue[j,] <- cov3$coverage
    DP.iobstrue[j] <- cov3$DP
    
    cov4 <- coverage(theta.true,c(est.lambda[,j],est.alpha[,j]),
                     Iobs.theta.est[,,j],rate,n)
    coverage.iobsest[j,] <- cov4$coverage
    DP.iobsest[j] <- cov4$DP
  }
}

### Computation of the empirical covering rates
colMeans(coverage.fisher.mcmc)
colMeans(coverage.iscotrue)
colMeans(coverage.iscoest)
colMeans(coverage.iobstrue[which(DP.iobstrue==1),])
colMeans(coverage.iobsest[which(DP.iobsest==1),])

## Graphical representation of the empirical distribution of the FIM estimates

DataRes <- data.frame(EstF11=c(sqrt(n)*(Isco.true.theta[1,1,]-fisher.mcmc[1,1]),
                               sqrt(n)*(Iobs.true.theta[1,1,]-fisher.mcmc[1,1])),
                      EstF22=c(sqrt(n)*(Isco.true.theta[2,2,]-fisher.mcmc[2,2]),
                               sqrt(n)*(Iobs.true.theta[2,2,]-fisher.mcmc[2,2])),
                      EstF33=c(sqrt(n)*(Isco.true.theta[3,3,]-fisher.mcmc[3,3]),
                               sqrt(n)*(Iobs.true.theta[3,3,]-fisher.mcmc[3,3])),
                      EstF44=c(sqrt(n)*(Isco.true.theta[4,4,]-fisher.mcmc[4,4]),
                               sqrt(n)*(Iobs.true.theta[4,4,]-fisher.mcmc[4,4])),
                      EstF55=c(sqrt(n)*(Isco.true.theta[5,5,]-fisher.mcmc[5,5]),
                               sqrt(n)*(Iobs.true.theta[5,5,]-fisher.mcmc[5,5])),
                      EstF12=c(sqrt(n)*(Isco.true.theta[1,2,]-fisher.mcmc[1,2]),
                               sqrt(n)*(Iobs.true.theta[1,2,]-fisher.mcmc[1,2])),
                      EstF13=c(sqrt(n)*(Isco.true.theta[1,3,]-fisher.mcmc[1,3]),
                               sqrt(n)*(Iobs.true.theta[1,3,]-fisher.mcmc[1,3])),
                      EstF23=c(sqrt(n)*(Isco.true.theta[2,3,]-fisher.mcmc[2,3]),
                               sqrt(n)*(Iobs.true.theta[2,3,]-fisher.mcmc[2,3])),
                      EstF35=c(sqrt(n)*(Isco.true.theta[3,5,]-fisher.mcmc[3,5]),
                               sqrt(n)*(Iobs.true.theta[3,5,]-fisher.mcmc[3,5])),
                      EstF34=c(sqrt(n)*(Isco.true.theta[3,4,]-fisher.mcmc[3,4]),
                               sqrt(n)*(Iobs.true.theta[3,4,]-fisher.mcmc[3,4])),
                      EstF25=c(sqrt(n)*(Isco.true.theta[2,5,]-fisher.mcmc[2,5]),
                               sqrt(n)*(Iobs.true.theta[2,5,]-fisher.mcmc[2,5])),
                      EstF24=c(sqrt(n)*(Isco.true.theta[2,4,]-fisher.mcmc[2,4]),
                               sqrt(n)*(Iobs.true.theta[2,4,]-fisher.mcmc[2,4])),
                      EstF15=c(sqrt(n)*(Isco.true.theta[1,5,]-fisher.mcmc[1,5]),
                               sqrt(n)*(Iobs.true.theta[1,5,]-fisher.mcmc[1,5])),
                      EstF14=c(sqrt(n)*(Isco.true.theta[1,4,]-fisher.mcmc[1,4]),
                               sqrt(n)*(Iobs.true.theta[1,4,]-fisher.mcmc[1,4])),
                      Estimate=c(rep('I N,sco',nbsim),rep('I N,obs',nbsim)))


F11 <- ggplot(DataRes, aes(EstF11, color=Estimate)) +
  geom_density(bw=0.2,alpha=0.6) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  xlim(-1,1) +
  ggtitle(bquote('('~lambda[1]~','~ lambda[1]~')')) +
  theme(legend.position = c(-0.05, 0.95), plot.title = element_text(size=20,face="bold"))

F22 <- ggplot(DataRes, aes(EstF22, color=Estimate)) +
  geom_density(bw=0.2,alpha=0.6) +
  scale_fill_manual(values = c("#984EA3",'#E69F00'))  +
  xlab("") +
  ylab("") +
  xlim(-0.75,0.75) +
  ggtitle(bquote('('~lambda[2]~','~ lambda[2]~')')) +
  theme(legend.position = "none", plot.title = element_text(size=20,face="bold"))

F25 <- ggplot(DataRes, aes(EstF25, color=Estimate)) +
  geom_density(bw=0.2,alpha=0.6) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  xlim(-1.5,1.5) +
  ggtitle(bquote('('~lambda[2]~','~ alpha[2]~')')) +
  theme(legend.position = "none", plot.title = element_text(size=20,face="bold"))


plot_grid(F11, F22, F25, ncol = 3, nrow = 1)
