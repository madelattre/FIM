## R script for the numerical study in the linear mixed effects model
## ------------------------------------------------------------------

nsim <- 500            # number of replicates
seq.n <- c(20,100,500) # number of individuals 
j    <- 12             # number of observations per individual

## parameter values
beta   <- 3
sigma2 <- 5 
eta2   <- 2
theta.true <- matrix(c(beta,eta2,sigma2),ncol=1)

## computation of the exact FIM

fisher <- Fisher_LMM(beta,sigma2,eta2,j)

## loop executing the nsim replicates of the experiment

for (n in seq.n){
  
  ## definition of the objects to store the estimation of the Fisher Information
  ## matrix 
  
  resIobs.theta.true <- array(NA,dim=c(3,3,nsim))
  resIsco.theta.true <- array(NA,dim=c(3,3,nsim))
  resIobs.theta.est  <- array(NA,dim=c(3,3,nsim))
  resIsco.theta.est  <- array(NA,dim=c(3,3,nsim))
  
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
  
  ## counters for the estimation of the coverage rates
  coverage.iobs.theta.est <- 0
  coverage.isco.theta.est <- 0
  coverage.iobs.theta.true <- 0
  coverage.isco.theta.true <- 0
  coverage.true.fisher <- 0
  
  for (k in 1:nsim){
    
    ## data simulation
    random    <- rnorm(n,0,sqrt(eta2))
    residual  <- rnorm(n*j,0,sqrt(sigma2))
    randompop <- rep(random,j)
    id        <- rep(seq(1,n),j)
    obs       <- beta+randompop+residual
    datamat   <- matrix(obs,n,j)
    
    ## computation of the FIM estimators when knowing the true parameter values --
    
    resIobs.theta.true[,,k] <- Iobs_LMM(datamat,beta,sigma2,eta2)
    resIsco.theta.true[,,k] <- Isco_LMM(datamat,beta,sigma2,eta2)
    
    ## computation of the FIM estimators together with the parameter values ----
    
    est.mle   <- lmer(obs~(1|id),REML = F)
    variances <- as.data.frame(VarCorr(est.mle))
    beta.est   <- c(beta.est,est.mle@beta)
    eta2.est   <- c(eta2.est,variances[1,'vcov'])
    sigma2.est <- c(sigma2.est,variances[1,'vcov'])
    
    resIobs.theta.est[,,k] <- Iobs_LMM(datamat,est.mle@beta,variances[2,'vcov'],
                                       variances[1,'vcov'])
    resIsco.theta.est[,,k] <- Isco_LMM(datamat,est.mle@beta,variances[2,'vcov'],
                                       variances[1,'vcov'])
    

    
  }
  
  EstF11.true <- c(EstF11.true,c(resIsco.theta.true[1,1,]-fisher[1,1],
                       resIobs.theta.true[1,1,]-fisher[1,1]))
  EstF22.true <- c(EstF22.true,c(resIsco.theta.true[2,2,]-fisher[2,2],
                       resIobs.theta.true[2,2,]-fisher[2,2]))
  EstF33.true <- c(EstF33.true,c(resIsco.theta.true[3,3,]-fisher[3,3],
                       resIobs.theta.true[3,3,]-fisher[3,3]))
  EstF12.true <- c(EstF12.true,c(resIsco.theta.true[1,2,]-fisher[1,2],
                       resIobs.theta.true[1,2,]-fisher[1,2]))
  EstF13.true <- c(EstF13.true,c(resIsco.theta.true[1,3,]-fisher[1,3],
                       resIobs.theta.true[1,3,]-fisher[1,3]))
  EstF23.true <- c(EstF23.true,c(resIsco.theta.true[2,3,]-fisher[2,3],
                       resIobs.theta.true[2,3,]-fisher[2,3]))
  
  
  EstF11.est <- c(EstF11.est,c(resIsco.theta.est[1,1,]-fisher[1,1],
                                 resIobs.theta.est[1,1,]-fisher[1,1]))
  EstF22.est <- c(EstF22.est,c(resIsco.theta.est[2,2,]-fisher[2,2],
                                 resIobs.theta.est[2,2,]-fisher[2,2]))
  EstF33.est <- c(EstF33.est,c(resIsco.theta.est[3,3,]-fisher[3,3],
                                 resIobs.theta.est[3,3,]-fisher[3,3]))
  EstF12.est <- c(EstF12.est,c(resIsco.theta.est[1,2,]-fisher[1,2],
                                 resIobs.theta.est[1,2,]-fisher[1,2]))
  EstF13.est <- c(EstF13.est,c(resIsco.theta.est[1,3,]-fisher[1,3],
                                 resIobs.theta.est[1,3,]-fisher[1,3]))
  EstF23.est <- c(EstF23.est,c(resIsco.theta.est[2,3,]-fisher[2,3],
                                 resIobs.theta.est[2,3,]-fisher[2,3]))
  
}

DataRes <- data.frame(EstF11.true=EstF11.true, EstF22.true=EstF22.true, 
                      EstF33.true=EstF33.true, EstF12.true=EstF12.true, 
                      EstF13.true=EstF13.true, EstF23.true=EstF23.true,
                      EstF11.est=EstF11.est, EstF22.est=EstF22.est, 
                      EstF33.est=EstF33.est, EstF12.est=EstF12.est, 
                      EstF13.est=EstF13.est, EstF23.est=EstF23.est,
                      beta.est=beta.est, eta2.est=eta2.est, sigma2.est=sigma2.est,
                      Estimate=rep(c(rep('I n,sco',nsim),rep('I n,obs',nsim)),
                                   length(seq.n)),
                      n=rep(seq.n,each=nsim*2)
)

save(DataRes,file="Rfiles/simusLMM.Rdata")



