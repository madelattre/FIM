load("Rfiles/ResNLMEnonexponential.Rdata")

vpop     <- 31
kapop    <- 1.6
clpop    <- 2.8
omega2ka <- 0.40
omega2cl <- 0.40
sigma2   <- 0.75

theta.true <- matrix(c(kapop,vpop,clpop,omega2ka,omega2cl,sigma2),6,1)


n     <- 100
rate  <- 0.95 

covering.isco <- rep(0,6)
covering.iobs <- rep(0,6)


nbsim <- dim(resNLMEnonExp$thetaest)[2]

for (k in 1:500){
  ICinfisco <- resNLMEnonExp$thetaest[,k] - 
    qnorm(1-(1-rate)/2,0,1)*sqrt(diag(solve(resNLMEnonExp$isco[,,k]))/n)
  ICsupisco <- resNLMEnonExp$thetaest[,k] + 
    qnorm(1-(1-rate)/2,0,1)*sqrt(diag(solve(resNLMEnonExp$isco[,,k]))/n)
  
  covering.isco <- covering.isco+as.numeric((theta.true<=ICsupisco)&(theta.true>=ICinfisco))
  
  ICinfiobs <- resNLMEnonExp$thetaest[,k] - 
    qnorm(1-(1-rate)/2,0,1)*sqrt(diag(solve(resNLMEnonExp$iobs[,,k]))/n)
  ICsupiobs <- resNLMEnonExp$thetaest[,k] + 
    qnorm(1-(1-rate)/2,0,1)*sqrt(diag(solve(resNLMEnonExp$iobs[,,k]))/n)
  
  covering.iobs <- covering.iobs+as.numeric((theta.true<=ICsupiobs)&(theta.true>=ICinfiobs))
}
