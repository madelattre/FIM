## Comparison with the numerical example of Meng & Spall (2017) on a mixture of two gaussian distributions
## The Fisher information matrix is estimated using the estimator based on the score, 
## with parameter values previously estimated via an EM algorithm
## (section 4.3 of the article)

rm(list=ls())

### EM algorithm

EM <- function(y, nbiter = 100){
  
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

### Computation of empirical covering rates and computation of the mean estimator based on the score based on a large number of simulated datasets


rate <- 0.95 # nominal rate


# True paramater values
probtrue <- 2/3 # mixture proportion
m1true   <- 3   # mean of the first mixture proportion
m2true   <- 0   # mean of the second mixture proportion


nrep <- 10000   # number of simulated datasets 

# Intermediary R objects
recouvprob <- 0
recouvm1   <- 0
recouvm2   <- 0
isco <- list()

## Start of the loop
## nrep data simulation + parameter and FIM estimation
for (j in 1:nrep){
  
  n <- 750
  
  z <- 1+rbinom(n,1,probtrue)
  y <- c()
  for (i in 1:n){
    y <- c(y,rnorm(1,m1true*(z[i]==1)+m2true*(z[i]==2)))
  }
  
  est <- EM(y)
  
  m1est   <- est$m1
  m2est   <- est$m2
  probest <- est$prob
  
  iscoest <- matrix(0,nrow=3,ncol=3)
  
  for (i in 1:n){
    denomi <- (1-probest)*exp(-1/2*(y[i]-m1est)^2) + probest*exp(-1/2*(y[i]-m2est)^2)
    espcondi <- as.matrix(1/denomi * c(exp(-1/2*(y[i]-m2est)^2)-exp(-1/2*(y[i]-m1est)^2),(y[i]-m1est)*(1-probest)*exp(-1/2*(y[i]-m1est)^2),(y[i]-m2est)*probest*exp(-1/2*(y[i]-m2est)^2)),nrow=3,ncol=1)
    iscoest <- iscoest + espcondi%*%t(espcondi)
  }

  ICinf <- c(probest,m1est,m2est) - qnorm(1-(1-rate)/2)*sqrt(diag(solve(iscoest)))
  ICsup <- c(probest,m1est,m2est) + qnorm(1-(1-rate)/2)*sqrt(diag(solve(iscoest)))
  
  if ((probtrue>=ICinf[1])&(probtrue<=ICsup[1])){recouvprob <- recouvprob + 1}
  if ((m1true>=ICinf[2])&(m1true<=ICsup[2])){recouvm1 <- recouvm1 + 1}
  if ((m2true>=ICinf[3])&(m2true<=ICsup[3])){recouvm2 <- recouvm2 + 1}
  
  isco[[j]] <- iscoest
}
## End of the loop

## Results

### empirical covering rates
recouvprob/nrep # mixture proportion
recouvm1/nrep   # mean of the first mixture component
recouvm2/nrep   # mean of the second mixture component

### empirical mean of the nrep estimates of the FIM

isco11 <- 0
isco12 <- 0
isco13 <- 0
isco22 <- 0
isco23 <- 0
isco33 <- 0

for (j in 1:nrep){
  isco11 <- isco11 + isco[[j]][1,1]
  isco12 <- isco12 + isco[[j]][1,2]
  isco13 <- isco13 + isco[[j]][1,3]
  isco22 <- isco22 + isco[[j]][2,2]
  isco23 <- isco23 + isco[[j]][2,3]
  isco33 <- isco33 + isco[[j]][3,3]
}

isco.mean <- matrix(c(isco11/nrep,isco12/nrep,isco13/nrep,isco12/nrep,isco22/nrep,isco23/nrep,isco13/nrep,isco23/nrep,isco33/nrep),ncol=3,nrow=3)
isco.mean

