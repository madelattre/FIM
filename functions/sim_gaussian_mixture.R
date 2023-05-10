## Simulation of observations according to a mixture of two Gaussian 
## distributions 

sim_gaussian_mixture <- function(n,m1,m2,prob){
  # n    : number of observations to be simulated
  # m1   : mean of the first Gaussian distribution
  # m2   : mean of the second Gaussian distribution
  # prob : mixture proportion of the second distribution
  
  z <- 1+rbinom(n,1,prob) # latent variable
  y <- c()                # observed variable
  
  for (i in 1:n){
    y <- c(y,rnorm(1,m1*(z[i]==1)+m2*(z[i]==2)))
  }
  
  return(y)
}

