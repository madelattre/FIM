## Monte-Carlo estimation of the true Fisher information matrix based on a large sample

source('fonctions/sim_poisson_mixture.R')
source('fonctions/sim_poisson_mixture.R')

nMC <- 100000000

alpha <- c(0.3,0.5) # mixture weights of the first K-1 components 
lambda <- c(2,5,9)  # parameter values of the K Poisson distribution of the mixture


y <- sim_poisson_mixture(nMC,lambda,alpha)
trueFIM <-  fisher_estimation_poisson_mixture(y, nMC, lambda, alpha)
trueFIM$Iobs
trueFIM$Isco
trueFIM <- (trueFIM$Isco+trueFIM$Iobs)/2
save(trueFIM,file='Rfiles/PoissonMixtureTrueFIM.Rdata')
