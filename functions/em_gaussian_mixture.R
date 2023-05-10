## EM algorithm for parameter estimation in the mixture of two Gaussian
## distributions with unknown means and known variances equal to 1

em_gaussian_mixture <- function(y, nbiter = 100){
  
  # y      : vector of observations
  # nbiter : number of em iterations
  
  # initialization using kmeans
  km   <- kmeans(y,centers=2)
  init <- list(Ppi = km$size[1]/sum(km$size), mu1 = km$centers[1], mu2 = km$centers[2])
  
  Ppi    <- rep(NA,(nbiter+1)) # estimations of the mixture proportion
  Ppi[1] <- init$Ppi
  mu1    <- rep(NA,(nbiter+1)) # mean estimations of the first mixture component
  mu1[1] <- init$mu1
  mu2    <- rep(NA,(nbiter+1)) # mean estimations of the second mixture component
  mu2[1] <- init$mu2
  
  for (k in 1: nbiter){ 
    
    ## E-step
    
    omega     <- matrix(0,ncol=2,nrow=n)
    omega[,1] <- dnorm(y,mean=mu1[k],sd=1)*(1-Ppi[k])/(dnorm(y,mean=mu1[k],sd=1)*(1-Ppi[k])+dnorm(y,mean=mu2[k],sd=1)*Ppi[k])
    omega[,2] <- dnorm(y,mean=mu2[k],sd=1)*Ppi[k]/(dnorm(y,mean=mu1[k],sd=1)*(1-Ppi[k])+dnorm(y,mean=mu2[k],sd=1)*Ppi[k])
    
    ## M-step
    mu1[k+1] <- sum(omega[,1]*y)/sum(omega[,1])
    mu2[k+1] <- sum(omega[,2]*y)/sum(omega[,2])
    Ppi[k+1] <- 1/n*sum(omega[,2])
  }
  
  if (mu1[nbiter+1]<mu2[nbiter+1]){
    m1 <- mu2[nbiter+1]
    m2 <- mu1[nbiter+1]
    prob <- 1-Ppi[nbiter+1]
  } else{
    m1 <- mu1[nbiter+1]
    m2 <- mu2[nbiter+1]
    prob <- Ppi[nbiter+1]
  }
  
  res = list(m1=m1,m2=m2,prob=prob)
  
  return(res)
}
