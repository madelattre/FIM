## Numerical study in the Gaussian mixture model from (Meng and Spall, 2017)
## This script computes the estimator Isco on a large number of simulated 
## datasets of size n=750. 


# True paramater values
probtrue <- 2/3 # mixture proportion
m1true   <- 3   # mean of the first mixture proportion
m2true   <- 0   # mean of the second mixture proportion

# Sample size
n <- 750 

# Number of simulated datasets
nrep <- 10000 

# Nominal rate for the computation of empirical coverage rates
rate <- 0.95 
 
# Intermediary R objects to store the results
recouvprob <- 0
recouvm1   <- 0
recouvm2   <- 0

isco <- array(NA,dim=c(nrep,3,3))

for (j in 1:nrep){
  
  y <- sim_gaussian_mixture(n,m1true,m2true,probtrue)
  
  est <- em_gaussian_mixture(y)
  
  iscoest <- fisher_estimation_gaussian_mixture(y,est$m1,est$m2,est$prob)
  
  ICinf <- 
    c(est$prob,est$m1,est$m2) - qnorm(1-(1-rate)/2)*sqrt(diag(solve(iscoest)))
  ICsup <- 
    c(est$prob,est$m1,est$m2) + qnorm(1-(1-rate)/2)*sqrt(diag(solve(iscoest)))
  
  if ((probtrue>=ICinf[1])&(probtrue<=ICsup[1])){recouvprob <- recouvprob + 1}
  if ((m1true>=ICinf[2])&(m1true<=ICsup[2])){recouvm1 <- recouvm1 + 1}
  if ((m2true>=ICinf[3])&(m2true<=ICsup[3])){recouvm2 <- recouvm2 + 1}
  
  isco[j,,] <- iscoest
}

res <- list(isco=round(apply(isco,c(2,3),mean),3),recouvprob=recouvprob/nrep,
            recouvm1=recouvm1/nrep,recouvm2=recouvm2/nrep)
save(res,file="Rfiles/ResGaussianMixture.Rdata")
