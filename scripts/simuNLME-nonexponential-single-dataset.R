## Nonlinear mixed-effects model not belonging to the curved exponential family
## Estimation and convergence graphs on one simulated dataset

## Sample characteristics 

n     <- 100                             # number of subjects
times <- c(0.25,0.5,1,2,3.5,5,7,9,12,24) # observation times
j     <- length(times)                   # number of observations per subject
dose  <- 320                             # dose

## True parameter values
vpop     <- 31
kapop    <- 1.6
clpop    <- 2.8
omega2ka <- 0.40
omega2cl <- 0.40 
sigma2   <- 0.75

## Simulation of individual parameters
vi  <- rep(vpop,n)
kai <- exp(rnorm(n,log(kapop),sd=sqrt(omega2ka)))
cli <- exp(rnorm(n,log(clpop),sd=sqrt(omega2cl)))

## Simulation of the observations
ypred <- c()
for (k in 1:n){
  ypred <- c(ypred,model1cptsim(cbind(kai,vi,cli),k,times,dose))
}

y <- ypred + rnorm(n*j,0,sd=sqrt(sigma2))

datasim <- data.frame(y=y,dose=rep(dose,n*j),time=rep(times,n),
                      subject=kronecker(1:n, rep(1,j)))

## Estimation

# Algorithmic settings
nbiterem     <- 3000
nbiterburnin <- 1500
theta0       <- list(vpop=(vpop-5)*runif(1,0.8,1.2),kapop=kapop*runif(1,0.8,1.2),
                     clpop=clpop*runif(1,0.8,1.2),omega2ka=omega2ka*runif(1,0.4,2),
                     omega2cl=omega2cl*runif(1,0.4,2), sigma2=sigma2*runif(1,0.4,2))

# Parameter and FIM estimation 
res <- saem_non_exp(datasim, nbiterem, nbiterburnin, theta0)

# Save one in 10 iterations due to limited space on Github
res$isco     <- res$isco[,,c(1,seq(10,nbiterem,10))]
res$thetaest <- res$thetaest[,c(1,seq(10,nbiterem,10))]

save(res,file="Rfiles/saem-non-exp-conv-plot.Rdata")

