### Numerical study in a linear mixed effects model : statistical properties of two estimators of the Fisher Information matrix 
### (section 4.1 of the article)

rm(list=ls())

## Source useful libraries ----------------------------------------------------

library(ggplot2)
library(cowplot)
#library(devtools)
#library(RColorBrewer)

## Defining R functions for computing the exact Fisher Information matrix, Iobs 
## and Isco in the linear mixed effects model ---------------------------------

## Here it is assumed that all individuals have the same number of observations.

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


## number of replicates
nsim <- 500

## sample size
# number of individuals
n <- 100
# number of observations per individual
j <- 12

## parameter values

beta <- 3
sigma2 <- 5 
eta2 <- 2

## computation of the exact FIM

fisher <- Fisher.LMM(beta,sigma2,eta2,j)


## Estimating the FIM when knowing the true parameter values ------------------

## definition of the objects to store the estimation of the Fisher Information matrix 

resIobs.theta.true <- array(NA,dim=c(3,3,nsim))
resIsco.theta.true <- array(NA,dim=c(3,3,nsim))
resIobs.theta.est <- array(NA,dim=c(3,3,nsim))
resIsco.theta.est <- array(NA,dim=c(3,3,nsim))




## loop executing the nsim replicates of the experiment

for (k in 1:nsim){
    
  ## data simulation
  random <- rnorm(n,0,sqrt(eta2))
  residual <- rnorm(n*j,0,sqrt(sigma2))
  randompop <- rep(random,j)
  id <- rep(seq(1,n),j)
  obs <- beta+randompop+residual
  datamat<-matrix(obs,n,j)
  
  
  
  ## computation of the FIM estimators when knowing the true parameter values --
  
  resIobs.theta.true[,,k] <- Iobs.LMM(datamat,beta,sigma2,eta2)
  resIsco.theta.true[,,k] <- Isco.LMM(datamat,beta,sigma2,eta2)
  
  ## computation of the FIM estimators together with the parameter values ----
  
  ## Estimation of the parameters
  est.mle   <- lmer(obs~(1|id),REML = F)
  variances <- as.data.frame(VarCorr(est.mle))
  
  ## Estimation of the FIM

  resIobs.theta.est[,,k] <- Iobs.LMM(datamat,est.mle@beta,variances[2,'vcov'],variances[1,'vcov'])
  resIsco.theta.est[,,k] <- Iobs.LMM(datamat,est.mle@beta,variances[2,'vcov'],variances[1,'vcov'])
  
  ## Coverages
  ## TBA
}

## TITLE TBA

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

