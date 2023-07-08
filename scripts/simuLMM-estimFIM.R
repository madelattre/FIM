## 1- R script for studying the asymptotic properties of Iobs and Isco in the 
## linear mixed effects model
## -----------------------------------------------------------------------------

library(lme4)

nsim  <- 500            # number of replicates
seq.n <- c(20,100,500)  # number of individuals 
j     <- 12             # number of observations per individual

## parameter values
beta   <- 3
sigma2 <- 5 
eta2   <- 2
theta.true <- matrix(c(beta,eta2,sigma2),ncol=1)

## R objects to store estimations

resIobs.theta.true <- array(NA,dim=c(3,3,nsim,length(seq.n)))
resIsco.theta.true <- array(NA,dim=c(3,3,nsim,length(seq.n)))
resIobs.theta.est  <- array(NA,dim=c(3,3,nsim,length(seq.n)))
resIsco.theta.est  <- array(NA,dim=c(3,3,nsim,length(seq.n)))

EstF11.true <- c()
EstF22.true <- c()
EstF33.true <- c()
EstF12.true <- c()
EstF13.true <- c()
EstF23.true <- c()

EstF11.est <- c()
EstF22.est <- c()
EstF33.est <- c()
EstF12.est <- c()
EstF13.est <- c()
EstF23.est <- c()

beta.est   <- c()
eta2.est   <- c()
sigma2.est <- c()

## loop executing the nsim replicates of the experiment

for (l in 1:length(seq.n)){
  
  n <- seq.n[l]
  
  beta.est.n   <- c()
  eta2.est.n   <- c()
  sigma2.est.n <- c()
  
  for (k in 1:nsim){
    
    ## data simulation
    random    <- rnorm(n,0,sqrt(eta2))
    residual  <- rnorm(n*j,0,sqrt(sigma2))
    randompop <- rep(random,j)
    id        <- rep(seq(1,n),j)
    obs       <- beta+randompop+residual
    datamat   <- matrix(obs,n,j)
    
    ## evaluation of the FIM estimators in the true parameter values 
    
    resIobs.theta.true[,,k,l] <- Iobs_LMM(datamat,beta,sigma2,eta2)
    resIsco.theta.true[,,k,l] <- Isco_LMM(datamat,beta,sigma2,eta2)
    
    ## evaluation of the FIM estimators in the estimated parameter values
    
    est.mle    <- lmer(obs~(1|id),REML = F)
    variances  <- as.data.frame(VarCorr(est.mle))
    beta.est.n   <- c(beta.est.n,est.mle@beta)
    eta2.est.n   <- c(eta2.est.n,variances[1,'vcov'])
    sigma2.est.n <- c(sigma2.est.n,variances[2,'vcov'])
    
    resIobs.theta.est[,,k,l] <- Iobs_LMM(datamat,est.mle@beta,
                                         variances[2,'vcov'],
                                       variances[1,'vcov'])
    resIsco.theta.est[,,k,l] <- Isco_LMM(datamat,est.mle@beta,
                                         variances[2,'vcov'],
                                       variances[1,'vcov'])
    
  }
  
  EstF11.true <- c(EstF11.true,c(resIsco.theta.true[1,1,,l],
                       resIobs.theta.true[1,1,,l]))
  EstF22.true <- c(EstF22.true,c(resIsco.theta.true[2,2,,l],
                       resIobs.theta.true[2,2,,l]))
  EstF33.true <- c(EstF33.true,c(resIsco.theta.true[3,3,,l],
                       resIobs.theta.true[3,3,,l]))
  EstF12.true <- c(EstF12.true,c(resIsco.theta.true[1,2,,l],
                       resIobs.theta.true[1,2,,l]))
  EstF13.true <- c(EstF13.true,c(resIsco.theta.true[1,3,,l],
                       resIobs.theta.true[1,3,,l]))
  EstF23.true <- c(EstF23.true,c(resIsco.theta.true[2,3,,l],
                       resIobs.theta.true[2,3,,l]))
  
  
  EstF11.est <- c(EstF11.est,c(resIsco.theta.est[1,1,,l],
                                 resIobs.theta.est[1,1,,l]))
  EstF22.est <- c(EstF22.est,c(resIsco.theta.est[2,2,,l],
                                 resIobs.theta.est[2,2,,l]))
  EstF33.est <- c(EstF33.est,c(resIsco.theta.est[3,3,,l],
                                 resIobs.theta.est[3,3,,l]))
  EstF12.est <- c(EstF12.est,c(resIsco.theta.est[1,2,,l],
                                 resIobs.theta.est[1,2,,l]))
  EstF13.est <- c(EstF13.est,c(resIsco.theta.est[1,3,,l],
                                 resIobs.theta.est[1,3,,l]))
  EstF23.est <- c(EstF23.est,c(resIsco.theta.est[2,3,,l],
                                 resIobs.theta.est[2,3,,l]))
  
  beta.est <- c(beta.est,rep(beta.est.n,2))
  eta2.est <- c(eta2.est,rep(eta2.est.n,2))
  sigma2.est <- c(sigma2.est,rep(sigma2.est.n,2))
  
}

DataRes <- data.frame(EstF11.true=EstF11.true, EstF22.true=EstF22.true, 
                      EstF33.true=EstF33.true, EstF12.true=EstF12.true, 
                      EstF13.true=EstF13.true, EstF23.true=EstF23.true,
                      EstF11.est=EstF11.est, EstF22.est=EstF22.est, 
                      EstF33.est=EstF33.est, EstF12.est=EstF12.est, 
                      EstF13.est=EstF13.est, EstF23.est=EstF23.est,
                      beta.est=beta.est, eta2.est=eta2.est, 
                      sigma2.est=sigma2.est,
                      Estimate=rep(c(rep('I n,sco',nsim),rep('I n,obs',nsim)),
                                   length(seq.n)),
                      n=rep(seq.n,each=nsim*2)
)

save(DataRes,file="Rfiles/simusLMM.Rdata")



