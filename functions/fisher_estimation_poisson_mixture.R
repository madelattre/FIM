## Function for computing Isco and Iobs for Fisher Information matrix estimation 
## in Poisson mixture models 

fisher_estimation_poisson_mixture <- function(y, lambda, alpha) {
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
    deriv1ind[k,]   <- 
      exp(-lambda[k])*lambda[k]^y*alpha[k]/denom*(y/lambda[k]-1)
    deriv1ind[K+k,] <- exp(-lambda[k])*lambda[k]^y/denom - 
      exp(-lambda[K])*lambda[K]^y/denom
  }
  deriv1ind[K,] <- 
    exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(y/lambda[K]-1)
  
  
  ## computation of conditional expectation of the second derivatives of the 
  ## complete data log-likelihood
  
  
  for (k in 1:(K-1)){
    deriv2[k,k]     <- 
      sum(exp(-lambda[k])*lambda[k]^y*alpha[k]/denom*(-y/lambda[k]^2))
    deriv2[K+k,K+k] <- 
      sum(-exp(-lambda[k])*lambda[k]^y/denom/alpha[k] -
                             exp(-lambda[K])*lambda[K]^y/denom*(1/(1-sum(alpha))))
  }
  
  for (k in 1:(K-2)){
    for (l in (k+1):(K-1)){
      deriv2[K+k,K+l] <- - 
        sum(exp(-lambda[K])*lambda[K]^y/denom*(1/(1-sum(alpha))))
      deriv2[K+l,K+k] <- deriv2[K+k,K+l] 
    }
  }
  
  deriv2[K,K]<-
    sum(exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(-y/lambda[K]^2))
  
  ## computation of the conditional covariance matrix of the first derivatives 
  ## of the complete data log-likelihood
  
  
  for (k in 1:(K-2)){
    covderiv[k,k] <- 
      sum(exp(-lambda[k])*lambda[k]^y*alpha[k]/denom*(-1+y/lambda[k])^2) 
    covderiv[k+K,k+K] <- sum(exp(-lambda[k])*lambda[k]^y/denom/alpha[k] +
                               exp(-lambda[K])*lambda[K]^y/denom/(1-sum(alpha))) 
    for (l in (k+1):(K-1)){
      covderiv[k+K,l+K] <- sum(exp(-lambda[K])*lambda[K]^y/denom/(1-sum(alpha))) 
      covderiv[l+K,k+K] <- covderiv[k+K,l+K]
    } 
    covderiv[k,K+k] <- sum(exp(-lambda[k])*lambda[k]^y/denom*(-1+y/lambda[k])) 
    covderiv[k+K,k] <- covderiv[k,K+k]
    
    covderiv[K,K+k] <- 
      sum(exp(-lambda[K])*lambda[K]^y/denom*(-1+y/lambda[K])*(-1)) 
    covderiv[K+k,K] <- covderiv[K,K+k]
  }
  
  covderiv[K-1,K-1] <- 
    sum(exp(-lambda[K-1])*lambda[K-1]^y*alpha[K-1]/denom*(-1+y/lambda[K-1])^2) 
  covderiv[2*K-1,2*K-1] <- sum(exp(-lambda[K-1])*lambda[K-1]^y/denom/alpha[K-1]+
                                 exp(-lambda[K])*lambda[K]^y/denom/(1-sum(alpha))) 
  covderiv[K-1,2*K-1] <- 
    sum(exp(-lambda[K-1])*lambda[K-1]^y/denom*(-1+y/lambda[K-1])) 
  covderiv[2*K-1,K-1] <- covderiv[K-1,2*K-1]
  
  covderiv[K,2*K-1] <- 
    sum(exp(-lambda[K])*lambda[K]^y/denom*(-1+y/lambda[K])*(-1)) 
  covderiv[2*K-1,K] <- covderiv[K,2*K-1]
  
  covderiv[K,K] <- 
    sum(exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(-1+y/lambda[K])^2) 
  
  ## computation of Isco and Iobs
  
  Isco <- deriv1ind%*%t(deriv1ind)/n
  Iobs <- deriv1ind%*%t(deriv1ind)/n - deriv2/n - covderiv/n
  
  
  res <- list(Isco = Isco, Iobs = Iobs)
  
  return(res)
}