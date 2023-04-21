## EM algorithm for parameter estimation in Poisson mixture models

em_poisson_mixture <- function(y,K,nbiter=100){
  # y : vector of observations
  # K : number of components of the mixture
  # nbiter : number of iterations of the algorithm
  
  n <- length(y)
  
  # initialization using kmeans
  
  km   <- kmeans(y,centers=K)
  init.alpha <- c(km$size[1]/sum(km$size), km$size[2]/sum(km$size), km$size[3]/sum(km$size))
  init.lambda <- c(km$centers[1], km$centers[2], km$centers[3])
  ord <- order(init.lambda)
  init.lambda <- init.lambda[ord]
  init.alpha <- init.alpha[ord][-K]
  
  lambda.est     <- matrix(NA,K,nbiter)
  alpha.est      <- matrix(NA,K-1,nbiter)
  lambda.est[,1] <- init.lambda
  alpha.est[,1]  <- init.alpha    
  
  prob <- matrix(NA,n,K-1) # posterior probabilities 
  
  for (l in 2:nbiter){
    
    denom <- 0
    for (k in 1:(K-1)){
      denom <- denom + 
        exp(-lambda.est[k,l-1])*lambda.est[k,l-1]^y*alpha.est[k,l-1]
    }
    denom <- denom + 
      exp(-lambda.est[K,l-1])*lambda.est[K,l-1]^y*(1-sum(alpha.est[,l-1]))
    
    ## E-step
    for (k in 1:(K-1)){
      prob[,k] <- exp(-lambda.est[k,l-1])*lambda.est[k,l-1]^y*alpha.est[k,l-1]/denom
    }
    
    ## M-step
    
    for (k in 1:(K-1)){
      lambda.est[k,l] <- sum(y*prob[,k])/sum(prob[,k])
      alpha.est[k,l] <- sum(prob[,k])/n
    }
    
    lambda.est[K,l] <- sum(y*(1-rowSums(prob)))/sum(1-rowSums(prob))
  }
  
  est <- list(lambda.est[,nbiter],alpha.est[,nbiter])
  
  return(est)
}