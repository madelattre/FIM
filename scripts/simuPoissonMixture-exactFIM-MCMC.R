## 2- R script for studying the asymptotic properties of Iobs and Isco in the 
## Poisson mixture model - Monte-Carlo estimation of the true Fisher information 
## matrix based on a very large sample 
## -----------------------------------------------------------------------------

nMC <- 100000000
alpha <- c(0.3,0.5) # mixture weights of the first K-1 components 
lambda <- c(2,5,9)  # parameter values of the K Poisson distributions 

y       <- sim_poisson_mixture(nMC,lambda,alpha)
trueFIM <-  fisher_estimation_poisson_mixture(y, nMC, lambda, alpha)
trueFIM <- (trueFIM$Isco+trueFIM$Iobs)/2
save(trueFIM,file='Rfiles/PoissonMixtureTrueFIM.Rdata')
