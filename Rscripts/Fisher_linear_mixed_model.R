### "Computing an empirical Fisher information matrix estimate in latent variable
### models through stochastic approximation."
### ---------------------------------------------------------------------------


### Numerical study in a linear mixed effects model : statistical properties of two estimators of the Fisher Information matrix 
### (section 4.1 of the article)

rm(list=ls())

## Source useful libraries ----------------------------------------------------

library(ggplot2)
library(cowplot)
library(lme4)

## Defining R functions for computing the exact Fisher Information matrix, Iobs 
## and Isco in the linear mixed effects model ---------------------------------

## It is assumed that all individuals have the same number of observations.

## Exact Fisher Information matrix

Fisher.LMM <- function(beta,sigma2,eta2,j){
  # beta    : value of the fixed-effect 
  # sigma2  : value of the residual variance 
  # eta2    : value of the random effects variance
  # j       : number of observations per individual
  
  crochet <- 2*j*eta2/(sigma2+eta2*j)^2/sigma2^3 +
              4/(sigma2+eta2*j)^2/sigma2^2 + 2/(sigma2+eta2*j)^3/sigma2
  alpha   <- j*(eta2+sigma2)/sigma2^3 - eta2/2*j*(j*eta2+sigma2)*crochet -
              (j-1)/2/sigma2^2 - 1/2/(sigma2+eta2*j)^2
  fisher  <- cbind(c(j/(sigma2+eta2*j),0,0),
                  c(0,j^2/2/(sigma2+eta2*j)^2,j/2/(sigma2+eta2*j)^2),
                  c(0,j/2/(sigma2+eta2*j)^2,alpha))

  return(fisher)
}

## Observed information matrix

Iobs.LMM <- function(datamat,beta,sigma2,eta2){
  # datamat : observations organized in a matrix. Each row is an individual
  #           vector of observations.
  # beta    : value of the fixed-effect 
  # sigma2  : value of the residual variance 
  # eta2    : value of the random effects variance
  
  n   <- dim(datamat)[1]
  j   <- dim(datamat)[2]
  obs <- as.vector(datamat)
  
  Iobs      <- matrix(0,3,3)
  Iobs[1,1] <- n*j/(sigma2+eta2*j)
  Iobs[2,1] <- j/(sigma2+eta2*j)^2*sum(obs-beta)
  Iobs[1,2] <- Iobs[2,1]
  Iobs[2,2] <- j/(sigma2+eta2*j)^3*sum(apply(datamat-beta,1,sum)^2)-n*j^2/2/(sigma2+eta2*j)^2
  Iobs[3,1] <- 1/(sigma2+eta2*j)^2*sum(obs-beta)
  Iobs[1,3] <- Iobs[3,1]
  Iobs[2,3] <- 1/(sigma2+eta2*j)^3*sum(apply(datamat-beta,1,sum)^2)-n*j/2/(sigma2+eta2*j)^2
  Iobs[3,2] <- Iobs[2,3]
  Iobs[3,3] <- 1/(sigma2)^3*sum((obs-beta)^2) -
                eta2*sum(apply(datamat-beta,1,sum)^2)*(j*eta2/(sigma2+eta2*j)^2/sigma2^3 +
                2/(sigma2+eta2*j)^2/sigma2^2 + 1/(sigma2+eta2*j)^3/sigma2 ) -
                n*(j-1)/2/sigma2^2 - n/2/(sigma2+eta2*j)^2
  Iobs <- Iobs/n
  
  return(Iobs)
}

## Estimator based on the score

Isco.LMM <- function(datamat,beta,sigma2,eta2){
  # datamat : observations organized in a matrix. Each row is an individual
  #           vector of observations.
  # beta    : value of the fixed-effect 
  # sigma2  : value of the residual variance 
  # eta2    : value of the random effects variance

  n   <- dim(datamat)[1]
  j   <- dim(datamat)[2]

  derivative     <- matrix(0,3,n)
  derivative[1,] <- apply(datamat-beta,1,sum)/(sigma2+eta2*j)
  derivative[2,] <- apply(datamat-beta,1,sum)^2/2/(sigma2+eta2*j)^2-j/2/(sigma2+eta2*j)
  derivative[3,] <- apply((datamat-beta)^2,1,sum)/2/sigma2^2 - 
    apply(datamat-beta,1,sum)^2 *eta2*(j*eta2+2*sigma2)/2/sigma2^2/(sigma2+eta2*j)^2 -
    1/2/(sigma2+eta2*j)-(j-1)/2/sigma2
  
  Isco <- derivative%*%t(derivative)/n
  
  return(Isco)
}

## Numerical experiments ------------------------------------------------------

nsim <- 1000 # number of replicates
n    <- 500 # number of individuals, either 20 or 100 or 500
j    <- 12 # number of observations per individual

## parameter values
beta   <- 3
sigma2 <- 5 
eta2   <- 2
theta.true <- matrix(c(beta,eta2,sigma2),ncol=1)

## computation of the exact FIM

fisher <- Fisher.LMM(beta,sigma2,eta2,j)

## definition of the objects to store the estimation of the Fisher Information
## matrix 

resIobs.theta.true <- array(NA,dim=c(3,3,nsim))
resIsco.theta.true <- array(NA,dim=c(3,3,nsim))
resIobs.theta.est <- array(NA,dim=c(3,3,nsim))
resIsco.theta.est <- array(NA,dim=c(3,3,nsim))

## counters for the estimation of the coverage rates
coverage.iobs.theta.est <- 0
coverage.isco.theta.est <- 0

coverage.iobs.theta.true <- 0
coverage.isco.theta.true <- 0

coverage.true.fisher <- 0


## loop executing the nsim replicates of the experiment

for (k in 1:nsim){
    
  ## data simulation
  random    <- rnorm(n,0,sqrt(eta2))
  residual  <- rnorm(n*j,0,sqrt(sigma2))
  randompop <- rep(random,j)
  id        <- rep(seq(1,n),j)
  obs       <- beta+randompop+residual
  datamat   <- matrix(obs,n,j)
   
  ## computation of the FIM estimators when knowing the true parameter values --
  
  resIobs.theta.true[,,k] <- Iobs.LMM(datamat,beta,sigma2,eta2)
  resIsco.theta.true[,,k] <- Isco.LMM(datamat,beta,sigma2,eta2)
  
  ## computation of the FIM estimators together with the parameter values ----
  
  ## Estimation of the parameters
  est.mle   <- lmer(obs~(1|id),REML = F)
  variances <- as.data.frame(VarCorr(est.mle))
  
  ## Estimation of the FIM

  resIobs.theta.est[,,k] <- Iobs.LMM(datamat,est.mle@beta,variances[2,'vcov'],
                                     variances[1,'vcov'])
  resIsco.theta.est[,,k] <- Isco.LMM(datamat,est.mle@beta,variances[2,'vcov'],
                                     variances[1,'vcov'])
  
  ## Coverage rates of confidence ellipsoids ...
  theta.est <- matrix(c(est.mle@beta,variances[1,'vcov'],variances[2,'vcov']),ncol=1)
  unit.vect <- matrix(rep(1,length(theta.est)),ncol=1)
  
  ## ... based on the true Fisher information matrix
  
  conf.inf.fisher.true <- theta.est -
    1.96/sqrt(n)*solve(chol(fisher))%*%unit.vect
  conf.sup.fisher.true <- theta.est +
    1.96/sqrt(n)*solve(chol(fisher))%*%unit.vect
  
  coverage.true.fisher <- coverage.true.fisher +
    (theta.true<=conf.sup.fisher.true)*(theta.true>=conf.inf.fisher.true)
  
  
  ## ... based on Isco computed in the MLE value of the parameters
  conf.inf.isco.theta.est <- theta.est -
    1.96/sqrt(n)*solve(chol(resIsco.theta.est[,,k]))%*%unit.vect
  conf.sup.isco.theta.est <- theta.est +
    1.96/sqrt(n)*solve(chol(resIsco.theta.est[,,k]))%*%unit.vect
  
  coverage.isco.theta.est <- coverage.isco.theta.est +
    (theta.true<=conf.sup.isco.theta.est)*(theta.true>=conf.inf.isco.theta.est)
  
  ## ... based on Isco computed in the true parameter values
  conf.inf.isco.theta.true <- theta.est -
    1.96/sqrt(n)*solve(chol(resIsco.theta.true[,,k]))%*%unit.vect
  conf.sup.isco.theta.true <- theta.est +
    1.96/sqrt(n)*solve(chol(resIsco.theta.true[,,k]))%*%unit.vect
  
  coverage.isco.theta.true <- coverage.isco.theta.true +
    (theta.true<=conf.sup.isco.theta.true)*(theta.true>=conf.inf.isco.theta.true)
  
  ## ... based on Iobs computed in the MLE value of the parameters
  conf.inf.iobs.theta.est <- theta.est -
    1.96/sqrt(n)*solve(chol(resIobs.theta.est[,,k]))%*%unit.vect
  conf.sup.iobs.theta.est <- theta.est +
    1.96/sqrt(n)*solve(chol(resIobs.theta.est[,,k]))%*%unit.vect
  
  coverage.iobs.theta.est <- coverage.iobs.theta.est +
    (theta.true<=conf.sup.iobs.theta.est)*(theta.true>=conf.inf.iobs.theta.est)
  
  ## ... based on Iobs computed in the true parameter values
  conf.inf.iobs.theta.true <- theta.est -
    1.96/sqrt(n)*solve(chol(resIobs.theta.true[,,k]))%*%unit.vect
  conf.sup.iobs.theta.true <- theta.est +
    1.96/sqrt(n)*solve(chol(resIobs.theta.true[,,k]))%*%unit.vect
  
  coverage.iobs.theta.true <- coverage.iobs.theta.true +
    (theta.true<=conf.sup.iobs.theta.true)*(theta.true>=conf.inf.iobs.theta.true)
  
}

## Print the coverage rates
coverage.isco.theta.true/nsim*100
coverage.isco.theta.est/nsim*100
coverage.iobs.theta.true/nsim*100
coverage.iobs.theta.est/nsim*100
coverage.true.fisher/nsim*100

## Graphical representation of the empirical densisities of the components of 
## the FIM estimators

DataRes <- data.frame(EstF11=c(sqrt(n)*(resIsco.theta.true[1,1,]-fisher[1,1]),
                               sqrt(n)*(resIobs.theta.true[1,1,]-fisher[1,1])),
                            EstF22=c(sqrt(n)*(resIsco.theta.true[2,2,]-fisher[2,2]),
                                     sqrt(n)*(resIobs.theta.true[2,2,]-fisher[2,2])),
                            EstF33=c(sqrt(n)*(resIsco.theta.true[3,3,]-fisher[3,3]),
                                     sqrt(n)*(resIobs.theta.true[3,3,]-fisher[3,3])),
                            EstF12=c(sqrt(n)*(resIsco.theta.true[1,2,]-fisher[1,2]),
                                     sqrt(n)*(resIobs.theta.true[1,2,]-fisher[1,2])),
                            EstF13=c(sqrt(n)*(resIsco.theta.true[1,3,]-fisher[1,3]),
                                     sqrt(n)*(resIobs.theta.true[1,3,]-fisher[1,3])),
                            EstF23=c(sqrt(n)*(resIsco.theta.true[2,3,]-fisher[2,3]),
                                     sqrt(n)*(resIobs.theta.true[2,3,]-fisher[2,3])),
                            Estimate=c(rep('I n,sco',nsim),rep('I n,obs',nsim)))

F22 <- ggplot(DataRes, aes(EstF22, fill=Estimate)) + 
  geom_density(bw=0.3,alpha=0.6) + 
  scale_fill_manual(values = c("#984EA3",'#E69F00')) + 
  xlab("") +
  ylab("") +
  xlim(-1.5,1.5) +
  ylim(0,2) +
  ggtitle(bquote('('~eta^2~','~eta^2~')')) +
  theme(legend.position = c(0,0.9), plot.title = element_text(size=20,face="bold"))

F33 <- ggplot(DataRes, aes(EstF33, fill=Estimate)) + 
  geom_density(bw=0.3,alpha=0.6) + 
  scale_fill_manual(values = c("#984EA3",'#E69F00')) + 
  xlab("") +
  ylab("") +
  xlim(-1.5,1.5) +
  ylim(0,2) +
  ggtitle(bquote('('~sigma^2~','~ sigma^2~')')) +
  theme(legend.position = "none", plot.title = element_text(size=20,face="bold"))


F12 <- ggplot(DataRes, aes(EstF12, fill=Estimate)) + 
  geom_density(bw=0.3,alpha=0.6) + 
  scale_fill_manual(values = c("#984EA3",'#E69F00')) + 
  xlab("") +
  ylab("") +
  xlim(-1.5,1.5) +
  ylim(0,2) +
  ggtitle(bquote('('~beta~','~ eta^2~')')) +
  theme(legend.position = "none", plot.title = element_text(size=20,face="bold"))

plot_grid(F22, F33, F12, ncol = 3, nrow = 1)

