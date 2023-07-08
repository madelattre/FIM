## 2- R script for studying the asymptotic properties of Iobs and Isco in the 
## Poisson mixture model
## -----------------------------------------------------------------------------

nbsim  <- 500           # number of replicates
alpha  <- c(0.3,0.5)    # mixture weights of the first K-1 components 
lambda <- c(2,5,9)      # parameter values of the K Poisson distributions
seq.n  <- c(20,100,500) # sample size

Iobs.theta.est  <- array(NA,dim=c(5,5,nbsim))
Isco.theta.est  <- array(NA,dim=c(5,5,nbsim))  
est.lambda      <- matrix(NA,3,nbsim)
est.alpha       <- matrix(NA,2,nbsim)

for (n in seq.n){
  for (j in 1:nbsim){
    
    ## Data simulation
    
    y <- sim_poisson_mixture(n,lambda,alpha)
    
    ## Parameter estimation
    
    em.est         <- em_poisson_mixture(y,3)
    est.lambda[,j] <- em.est[[1]]
    est.alpha[,j]  <- em.est[[2]]
    
    ## Computation of Isco and Iobs in the MLE value of the parameter
    
    res.theta.est       <- fisher_estimation_poisson_mixture(y, est.lambda[,j], 
                                                             est.alpha[,j])
    Iobs.theta.est[,,j] <- res.theta.est$Iobs
    Isco.theta.est[,,j] <- res.theta.est$Isco
  }
  
  ResSim   <- list(n=n,Isco=Isco.theta.est,Iobs=Iobs.theta.est,lambda=lambda,
                 alpha=alpha)
  filename <- paste('Rfiles/simusMixt_n',n,'.Rdata',sep="")
  save(ResSim,file=filename)
}






