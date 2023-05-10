## R model functions for the nonlinear PK model

## Functions returning the conditional mean of the observations given 
## the parameters psi and structural covariates, here dose and time

model1cptsim <- function(psi,id,tim,dose) {  
  # psi  : matrix of individual parameters (as many rows as individuals, 
  # as many columns as parameters)
  # id   : individual's number 
  # tim  : vector of observation times
  # dose : dose administered to the individual
  ka    <- psi[id,1]
  V     <- psi[id,2]
  CL    <- psi[id,3]
  k     <- CL/V
  ypred <- dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}


model1cpt<-function(psi,id,xidep) { 
  # psi   : matrix of individual parameters (as many rows as individuals, 
  # as many columns as parameters)
  # id    : individual's number 
  # xidep : matrix of doses (1st column) and observation times (2nd column)
  dose  <- xidep[,1]
  tim   <- xidep[,2]  
  ka    <- psi[id,1]
  V     <- psi[id,2]
  CL    <- psi[id,3]
  k     <- CL/V
  ypred <- dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

## Model function when the volume is fixed

model1cptV<-function(psi,id,xidep,V) {
  # psi  : matrix of individual parameters (as many rows as individuals, 
  # two columns)
  # id   : individual's number 
  # xidep: matrix of doses (1st column) and observation times (2nd column)
  # V    : value of the volume (fixed effect)
  dose  <- xidep[,1]
  tim   <- xidep[,2]
  ka    <- psi[id,1]
  CL    <- psi[id,2]
  k     <- CL/V
  ypred <- dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

## Derivative of the model function with respect to V

dVmodel1cpt<-function(psi,id,xidep,V) {
  # psi  : matrix of individual parameters (as many rows as individuals, 
  # two columns)
  # id   : individual's number 
  # xidep: matrix of doses (1st column) and observation times (2nd column)
  # V    : value of the volume (fixed effect)
  dose  <- xidep[,1]
  tim   <- xidep[,2]
  ka    <- psi[id,1]
  CL    <- psi[id,2]
  ypred <- -dose*ka^2/(V*ka-CL)^2*(exp(-CL/V*tim)-exp(-ka*tim))+
    dose*ka/(V*ka-CL)*CL/V^2*tim*exp(-CL/V*tim)
  return(ypred)
}

## Second derivative of the model function with respect to V

d2Vmodel1cpt<-function(psi,id,xidep,V) {
  # psi  : matrix of individual parameters (as many rows as individuals, 
  # two columns)
  # id   : individual's number 
  # xidep: matrix of doses (1st column) and observation times (2nd column)
  # V    : value of the volume (fixed effect)
  dose  <- xidep[,1]
  tim   <- xidep[,2]
  ka    <- psi[id,1]
  CL    <- psi[id,2]
  ypred <- 2*dose*ka^3/(V*ka-CL)^3*(exp(-CL/V*tim)-exp(-ka*tim))+
    dose*ka/(V*ka-CL)*CL/V^2*tim*exp(-CL/V*tim)*(-2*ka/(V*ka-CL)-2/V+CL*tim/V^2)
  return(ypred)
}