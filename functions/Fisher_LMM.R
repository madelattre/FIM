## Computation of the exact Fisher Information matrix in the linear mixed-effects
## model

Fisher_LMM <- function(beta,sigma2,eta2,j){
  # beta    : value of the fixed-effect 
  # sigma2  : value of the residual variance 
  # eta2    : value of the random effects variance
  # j       : number of observations per individual
  
  crochet <- 2*j*eta2/(sigma2+eta2*j)^2/sigma2^3 + 
    4/(sigma2+eta2*j)^2/sigma2^2 + 
    2/(sigma2+eta2*j)^3/sigma2
  
  alpha   <- j*(eta2+sigma2)/sigma2^3 - 
    eta2/2*j*(j*eta2+sigma2)*crochet -
    (j-1)/2/sigma2^2 - 
    1/2/(sigma2+eta2*j)^2
  
  fisher  <- cbind(c(j/(sigma2+eta2*j),0,0),
                   c(0,j^2/2/(sigma2+eta2*j)^2,j/2/(sigma2+eta2*j)^2),
                   c(0,j/2/(sigma2+eta2*j)^2,alpha))
  
  return(fisher)
}