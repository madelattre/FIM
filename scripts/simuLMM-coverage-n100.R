load("Rfiles/simusLMM.Rdata")

beta   <- 3
sigma2 <- 5 
eta2   <- 2
theta.true <- c(beta,eta2,sigma2)

DataRes.n500 <- DataRes %>% filter(n==100)
DataRes.n500.Isco <- DataRes.n500 %>% filter(Estimate=="I n,sco")
DataRes.n500.Iobs <- DataRes.n500 %>% filter(Estimate=="I n,obs")

seq.rate <- c(0.01,0.05,0.10)

nsim <- length(DataRes.n500.Iobs$beta.est)
n <- 100

coverage.iobs.theta.est <- matrix(0,length(seq.rate),3)
coverage.isco.theta.est <- matrix(0,length(seq.rate),3)
coverage.iobs.theta.true <- matrix(0,length(seq.rate),3)
coverage.isco.theta.true <- matrix(0,length(seq.rate),3)

for (rate in seq.rate){
  

  
  l <- which(rate==seq.rate)
  
  for (k in 1:nsim){
    
    theta.est <- c(DataRes.n500.Iobs$beta.est[k],DataRes.n500.Iobs$eta2.est[k],
                   DataRes.n500.Iobs$sigma2.est[k])
    
    Iobs.true <- matrix(c(DataRes.n500.Iobs$EstF11.true[k],DataRes.n500.Iobs$EstF12.true[k],
                          DataRes.n500.Iobs$EstF13.true[k],DataRes.n500.Iobs$EstF12.true[k],
                          DataRes.n500.Iobs$EstF22.true[k],DataRes.n500.Iobs$EstF23.true[k],
                          DataRes.n500.Iobs$EstF13.true[k],DataRes.n500.Iobs$EstF23.true[k],
                          DataRes.n500.Iobs$EstF33.true[k]),nrow=3)
    
    Isco.true <- matrix(c(DataRes.n500.Isco$EstF11.true[k],DataRes.n500.Isco$EstF12.true[k],
                          DataRes.n500.Isco$EstF13.true[k],DataRes.n500.Isco$EstF12.true[k],
                          DataRes.n500.Isco$EstF22.true[k],DataRes.n500.Isco$EstF23.true[k],
                          DataRes.n500.Isco$EstF13.true[k],DataRes.n500.Isco$EstF23.true[k],
                          DataRes.n500.Isco$EstF33.true[k]),nrow=3)
    
    Iobs.est <- matrix(c(DataRes.n500.Iobs$EstF11.est[k],DataRes.n500.Iobs$EstF12.est[k],
                         DataRes.n500.Iobs$EstF13.est[k],DataRes.n500.Iobs$EstF12.est[k],
                         DataRes.n500.Iobs$EstF22.est[k],DataRes.n500.Iobs$EstF23.est[k],
                         DataRes.n500.Iobs$EstF13.est[k],DataRes.n500.Iobs$EstF23.est[k],
                         DataRes.n500.Iobs$EstF33.est[k]),nrow=3)
    
    Isco.est <- matrix(c(DataRes.n500.Isco$EstF11.est[k],DataRes.n500.Isco$EstF12.est[k],
                         DataRes.n500.Isco$EstF13.est[k],DataRes.n500.Isco$EstF12.est[k],
                         DataRes.n500.Isco$EstF22.est[k],DataRes.n500.Isco$EstF23.est[k],
                         DataRes.n500.Isco$EstF13.est[k],DataRes.n500.Isco$EstF23.est[k],
                         DataRes.n500.Isco$EstF33.est[k]),nrow=3)
    
    IC.Iobs.true.inf <- (theta.est-1/sqrt(n)*qnorm(1-rate/2)*sqrt(diag(solve(Iobs.true))))
    IC.Iobs.true.sup <- (theta.est+1/sqrt(n)*qnorm(1-rate/2)*sqrt(diag(solve(Iobs.true))))
    
    IC.Isco.true.inf <- (theta.est-1/sqrt(n)*qnorm(1-rate/2)*sqrt(diag(solve(Isco.true))))
    IC.Isco.true.sup <- (theta.est+1/sqrt(n)*qnorm(1-rate/2)*sqrt(diag(solve(Isco.true))))
    
    IC.Iobs.est.inf <- (theta.est-1/sqrt(n)*qnorm(1-rate/2)*sqrt(diag(solve(Iobs.est))))
    IC.Iobs.est.sup <- (theta.est+1/sqrt(n)*qnorm(1-rate/2)*sqrt(diag(solve(Iobs.est))))
    
    IC.Isco.est.inf <- (theta.est-1/sqrt(n)*qnorm(1-rate/2)*sqrt(diag(solve(Isco.est))))
    IC.Isco.est.sup <- (theta.est+1/sqrt(n)*qnorm(1-rate/2)*sqrt(diag(solve(Isco.est))))
    
    coverage.iobs.theta.true[l,] <- coverage.iobs.theta.true[l,] + (theta.true>=IC.Iobs.true.inf)*(theta.true<=IC.Iobs.true.sup)
    coverage.isco.theta.true[l,] <- coverage.isco.theta.true[l,] + (theta.true>=IC.Isco.true.inf)*(theta.true<=IC.Isco.true.sup)
    coverage.iobs.theta.est[l,] <- coverage.iobs.theta.est[l,] + (theta.true>=IC.Iobs.est.inf)*(theta.true<=IC.Iobs.est.sup)
    coverage.isco.theta.est[l,] <- coverage.isco.theta.est[l,] + (theta.true>=IC.Isco.est.inf)*(theta.true<=IC.Isco.est.sup)
    
  }
  
}

dataCoverageLMM <- cbind(rate = rep(1-seq.rate,each=4),
                         fisher = rep(c('Isco','Isco','Iobs','Iobs'),3),
                         theta = rep(c('Known','Estimated','Known','Estimated'),3),
                         beta = c(coverage.isco.theta.true[1,1]/nsim,
                                  coverage.isco.theta.est[1,1]/nsim,
                                  coverage.iobs.theta.true[1,1]/nsim,
                                  coverage.iobs.theta.est[1,1]/nsim,
                                  coverage.isco.theta.true[2,1]/nsim,
                                  coverage.isco.theta.est[2,1]/nsim,
                                  coverage.iobs.theta.true[2,1]/nsim,
                                  coverage.iobs.theta.est[2,1]/nsim,
                                  coverage.isco.theta.true[3,1]/nsim,
                                  coverage.isco.theta.est[3,1]/nsim,
                                  coverage.iobs.theta.true[3,1]/nsim,
                                  coverage.iobs.theta.est[3,1]/nsim),
                         eta2 = c(coverage.isco.theta.true[1,2]/nsim,
                                  coverage.isco.theta.est[1,2]/nsim,
                                  coverage.iobs.theta.true[1,2]/nsim,
                                  coverage.iobs.theta.est[1,2]/nsim,
                                  coverage.isco.theta.true[2,2]/nsim,
                                  coverage.isco.theta.est[2,2]/nsim,
                                  coverage.iobs.theta.true[2,2]/nsim,
                                  coverage.iobs.theta.est[2,2]/nsim,
                                  coverage.isco.theta.true[3,2]/nsim,
                                  coverage.isco.theta.est[3,2]/nsim,
                                  coverage.iobs.theta.true[3,2]/nsim,
                                  coverage.iobs.theta.est[3,2]/nsim),
                         sigma2 = c(coverage.isco.theta.true[1,3]/nsim,
                                  coverage.isco.theta.est[1,3]/nsim,
                                  coverage.iobs.theta.true[1,3]/nsim,
                                  coverage.iobs.theta.est[1,3]/nsim,
                                  coverage.isco.theta.true[2,3]/nsim,
                                  coverage.isco.theta.est[2,3]/nsim,
                                  coverage.iobs.theta.true[2,3]/nsim,
                                  coverage.iobs.theta.est[2,3]/nsim,
                                  coverage.isco.theta.true[3,3]/nsim,
                                  coverage.isco.theta.est[3,3]/nsim,
                                  coverage.iobs.theta.true[3,3]/nsim,
                                  coverage.iobs.theta.est[3,3]/nsim))

dataCoverageLMM <- as.data.frame(dataCoverageLMM)
dataCoverageLMM <- flextable(dataCoverageLMM, cwidth=1.2)
dataCoverageLMM <- theme_box(dataCoverageLMM)
dataCoverageLMM <- bold(dataCoverageLMM,j=1)
dataCoverageLMM <- set_header_labels(
  x = dataCoverageLMM, values = c(rate="1-α", fisher="Fisher est.", theta="θ",
                                  beta="β", eta2="η2", sigma2="σ2"))
dataCoverageLMM <- dataCoverageLMM %>% merge_at(i=1:4,j=1) %>% merge_at(i=5:8,j=1) %>%
 merge_at(i=9:12,j=1) %>% merge_at(i=1:2,j=2) %>% merge_at(i=3:4,j=2) %>%
 merge_at(i=5:6,j=2) %>% merge_at(i=7:8,j=2) %>% merge_at(i=9:10,j=2) %>%
 merge_at(i=11:12,j=2)
dataCoverageLMM <- width(dataCoverageLMM, width = 0.8)
dataCoverageLMM