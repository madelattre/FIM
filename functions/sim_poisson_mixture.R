## Simulation of observations according to a Poisson mixture model

sim_poisson_mixture <- function(n,lambda,alpha){
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