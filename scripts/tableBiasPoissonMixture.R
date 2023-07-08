load('Rfiles/PoissonMixtureTrueFIM.Rdata')
load("Rfiles/simusMixt_n20.Rdata") 
biasIscoMixt20 <- round(apply(ResSim$Isco,c(1,2),mean) - trueFIM,5)
sdIscoMixt20 <- round(apply(ResSim$Isco,c(1,2),sd),5)
biasIobsMixt20 <- round(apply(ResSim$Iobs,c(1,2),mean) - trueFIM,5)
sdIobsMixt20 <- round(apply(ResSim$Iobs,c(1,2),sd),5)


load("Rfiles/simusMixt_n100.Rdata")
biasIscoMixt100 <- round(apply(ResSim$Isco,c(1,2),mean) - trueFIM,5)
sdIscoMixt100 <- round(apply(ResSim$Isco,c(1,2),sd),5)
biasIobsMixt100 <- round(apply(ResSim$Iobs,c(1,2),mean) - trueFIM,5)
sdIobsMixt100 <- round(apply(ResSim$Iobs,c(1,2),sd),5)

load("Rfiles/simusMixt_n500.Rdata")
biasIscoMixt500 <- round(apply(ResSim$Isco,c(1,2),mean) - trueFIM,5)
sdIscoMixt500 <- round(apply(ResSim$Isco,c(1,2),sd),5)
biasIobsMixt500 <- round(apply(ResSim$Iobs,c(1,2),mean) - trueFIM,5)
sdIobsMixt500 <- round(apply(ResSim$Iobs,c(1,2),sd),5)

dataBiasMixt <- cbind(component = c(1,2,3,4,5,6),
                      mIsco20 = c(diag(biasIscoMixt20)[2:5],biasIscoMixt20[2,3],biasIscoMixt20[3,5]),
                      mIobs20 = c(diag(biasIobsMixt20)[2:5],biasIobsMixt20[2,3],biasIobsMixt20[3,5]),
                      mIsco100 =c(diag(biasIscoMixt100)[2:5],biasIscoMixt100[2,3],biasIscoMixt100[3,5]),
                      mIobs100 =c(diag(biasIobsMixt100)[2:5],biasIobsMixt100[2,3],biasIobsMixt100[3,5]),
                      mIsco500 =c(diag(biasIscoMixt500)[2:5],biasIscoMixt500[2,3],biasIscoMixt500[3,5]),
                      mIobs500 =c(diag(biasIobsMixt500)[2:5],biasIobsMixt500[2,3],biasIobsMixt500[3,5]))


dataBiasMixt <- as.data.frame(dataBiasMixt)
dataBiasMixt$component <-factor(dataBiasMixt$component, levels=c(1,2,3,4,5,6),
                                labels=c("1"='(λ2,λ2)',
                                         "2"='(λ3,λ3)',
                                         "3"='(α1,α1)',
                                         "4"='(α2,α2)',
                                         "5"='(λ2,λ3)',
                                         "6"='(λ3,α2)'))
dataBiasMixt <- flextable(dataBiasMixt, cwidth=1.2)
dataBiasMixt <- add_header_row(
  x = dataBiasMixt, values = c("", "n=20", "n=100", "n=500"),
  colwidths = c(1, 2, 2, 2))
dataBiasMixt <- set_header_labels(dataBiasMixt, mIsco20="Isco", mIobs20="Iobs", 
                                  mIsco100="Isco", mIobs100="Iobs", mIsco500="Isco", 
                                  mIobs500="Iobs", component="")
dataBiasMixt <- align(dataBiasMixt, part = "all", align = "center")
dataBiasMixt <- bold(dataBiasMixt, j=c(2,4,6))
dataBiasMixt <- width(dataBiasMixt, width = 0.8)
dataBiasMixt