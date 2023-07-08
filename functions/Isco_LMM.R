## Computation of the estimator of the Fisher information matrix based on 
## the score

Isco_LMM <- function(datamat,beta,sigma2,eta2){
  # datamat : observations organized in a matrix. Each row is an individual
  #           vector of observations.
  # beta    : value of the fixed-effect 
  # sigma2  : value of the residual variance 
  # eta2    : value of the random effects variance
  
  n   <- dim(datamat)[1]
  j   <- dim(datamat)[2]
  
  derivative     <- matrix(0,3,n)
  derivative[1,] <- apply(datamat-beta,1,sum)/(sigma2+eta2*j)
  derivative[2,] <- apply(datamat-beta,1,sum)^2/2/(sigma2+eta2*j)^2-
    j/2/(sigma2+eta2*j)
  derivative[3,] <- apply((datamat-beta)^2,1,sum)/2/sigma2^2 - 
    apply(datamat-beta,1,sum)^2*eta2*(j*eta2+2*sigma2)/2/sigma2^2/(sigma2+eta2*j)^2 -
    1/2/(sigma2+eta2*j)-(j-1)/2/sigma2
  
  Isco <- derivative%*%t(derivative)/n
  
  return(Isco)
}