## Computation of the observed information matrix in the linear mixed-effects 
## model

Iobs_LMM <- function(datamat,beta,sigma2,eta2){
  # datamat : observations organized in a matrix. Each row is an individual
  #           vector of observations.
  # beta    : value of the fixed-effect 
  # sigma2  : value of the residual variance 
  # eta2    : value of the random effects variance
  
  n   <- dim(datamat)[1]
  j   <- dim(datamat)[2]
  obs <- as.vector(datamat)
  
  Iobs      <- matrix(0,3,3)
  Iobs[1,1] <- n*j/(sigma2+eta2*j)
  Iobs[2,1] <- j/(sigma2+eta2*j)^2*sum(obs-beta)
  Iobs[1,2] <- Iobs[2,1]
  Iobs[2,2] <- j/(sigma2+eta2*j)^3*sum(apply(datamat-beta,1,sum)^2)-
    n*j^2/2/(sigma2+eta2*j)^2
  Iobs[3,1] <- 1/(sigma2+eta2*j)^2*sum(obs-beta)
  Iobs[1,3] <- Iobs[3,1]
  Iobs[2,3] <- 1/(sigma2+eta2*j)^3*sum(apply(datamat-beta,1,sum)^2)-
    n*j/2/(sigma2+eta2*j)^2
  Iobs[3,2] <- Iobs[2,3]
  Iobs[3,3] <- 1/(sigma2)^3*sum((obs-beta)^2) -
    eta2*sum(apply(datamat-beta,1,sum)^2)*(j*eta2/(sigma2+eta2*j)^2/sigma2^3 +
    2/(sigma2+eta2*j)^2/sigma2^2 + 1/(sigma2+eta2*j)^3/sigma2 ) -
    n*(j-1)/2/sigma2^2 - n/2/(sigma2+eta2*j)^2
  Iobs <- Iobs/n
  
  return(Iobs)
}