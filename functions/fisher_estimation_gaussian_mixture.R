## Function for computing Isco for Fisher Information matrix estimation 
## in the mixture of two Gaussian distributions of variances 1

fisher_estimation_gaussian_mixture <- function(y, m1, m2, prob) {
  # y      : vector of observations
  # m1   : mean of the first Gaussian distribution
  # m2   : mean of the second Gaussian distribution
  # prob : mixture proportion of the second distribution
  
  isco <- matrix(0,nrow=3,ncol=3)
  
  n <- length(y)
  
  for (i in 1:n){
    denomi <- (1-prob)*exp(-1/2*(y[i]-m1)^2) + prob*exp(-1/2*(y[i]-m2)^2)
    espcondi <- as.matrix(1/denomi * c(exp(-1/2*(y[i]-m2)^2)-exp(-1/2*(y[i]-m1)^2),
                                       (y[i]-m1)*(1-prob)*exp(-1/2*(y[i]-m1)^2),
                                       (y[i]-m2)*prob*exp(-1/2*(y[i]-m2)^2)),
                          nrow=3,ncol=1)
    isco <- isco + espcondi%*%t(espcondi)
  }
  
  return(isco)
  
}