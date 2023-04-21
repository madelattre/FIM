load("Rfiles/simusLMM.Rdata")


DataRes.n20 <- DataRes %>% filter(n==20)
DataRes.n20.Isco <- DataRes.n20 %>% filter(Estimate=="I n,sco")
DataRes.n20.Iobs <- DataRes.n20 %>% filter(Estimate=="I n,obs")
mIsco20 <- round(c(mean(sqrt(20)*DataRes.n20.Isco$EstF11.true),mean(sqrt(20)*DataRes.n20.Isco$EstF22.true),
                   mean(sqrt(20)*DataRes.n20.Isco$EstF33.true),mean(sqrt(20)*DataRes.n20.Isco$EstF12.true),
                   mean(sqrt(20)*DataRes.n20.Isco$EstF13.true),mean(sqrt(20)*DataRes.n20.Isco$EstF23.true)),5)
mIobs20 <- round(c(mean(sqrt(20)*DataRes.n20.Iobs$EstF11.true),mean(sqrt(20)*DataRes.n20.Iobs$EstF22.true),
                   mean(sqrt(20)*DataRes.n20.Iobs$EstF33.true),mean(sqrt(20)*DataRes.n20.Iobs$EstF12.true),
                   mean(sqrt(20)*DataRes.n20.Iobs$EstF13.true),mean(sqrt(20)*DataRes.n20.Iobs$EstF23.true)),5)

sIsco20 <- round(c(sd(sqrt(20)*DataRes.n20.Isco$EstF11.true),sd(sqrt(20)*DataRes.n20.Isco$EstF22.true),
                   sd(sqrt(20)*DataRes.n20.Isco$EstF33.true),sd(sqrt(20)*DataRes.n20.Isco$EstF12.true),
                   sd(sqrt(20)*DataRes.n20.Isco$EstF13.true),sd(sqrt(20)*DataRes.n20.Isco$EstF23.true)),5)
sIobs20 <- round(c(sd(sqrt(20)*DataRes.n20.Iobs$EstF11.true),sd(sqrt(20)*DataRes.n20.Iobs$EstF22.true),
                   sd(sqrt(20)*DataRes.n20.Iobs$EstF33.true),sd(sqrt(20)*DataRes.n20.Iobs$EstF12.true),
                   sd(sqrt(20)*DataRes.n20.Iobs$EstF13.true),sd(sqrt(20)*DataRes.n20.Iobs$EstF23.true)),5)


DataRes.n100 <- DataRes %>% filter(n==100)
DataRes.n100.Isco <- DataRes.n100 %>% filter(Estimate=="I n,sco")
DataRes.n100.Iobs <- DataRes.n100 %>% filter(Estimate=="I n,obs")
mIsco100 <- round(c(mean(sqrt(100)*DataRes.n100.Isco$EstF11.true),mean(sqrt(100)*DataRes.n100.Isco$EstF22.true),
                    mean(sqrt(100)*DataRes.n100.Isco$EstF33.true),mean(sqrt(100)*DataRes.n100.Isco$EstF12.true),
                    mean(sqrt(100)*DataRes.n100.Isco$EstF13.true),mean(sqrt(100)*DataRes.n100.Isco$EstF23.true)),5)
mIobs100 <- round(c(mean(sqrt(100)*DataRes.n100.Iobs$EstF11.true),mean(sqrt(100)*DataRes.n100.Iobs$EstF22.true),
                    mean(sqrt(100)*DataRes.n100.Iobs$EstF33.true),mean(sqrt(100)*DataRes.n100.Iobs$EstF12.true),
                    mean(sqrt(100)*DataRes.n100.Iobs$EstF13.true),mean(sqrt(100)*DataRes.n100.Iobs$EstF23.true)),5)

sIsco100 <- round(c(sd(sqrt(100)*DataRes.n100.Isco$EstF11.true),sd(sqrt(100)*DataRes.n100.Isco$EstF22.true),
                    sd(sqrt(100)*DataRes.n100.Isco$EstF33.true),sd(sqrt(100)*DataRes.n100.Isco$EstF12.true),
                    sd(sqrt(100)*DataRes.n100.Isco$EstF13.true),sd(sqrt(100)*DataRes.n100.Isco$EstF23.true)),5)
sIobs100 <- round(c(sd(sqrt(100)*DataRes.n100.Iobs$EstF11.true),sd(sqrt(100)*DataRes.n100.Iobs$EstF22.true),
                    sd(sqrt(100)*DataRes.n100.Iobs$EstF33.true),sd(sqrt(100)*DataRes.n100.Iobs$EstF12.true),
                    sd(sqrt(100)*DataRes.n100.Iobs$EstF13.true),sd(sqrt(100)*DataRes.n100.Iobs$EstF23.true)),5)

DataRes.n500 <- DataRes %>% filter(n==500)
DataRes.n500.Isco <- DataRes.n500 %>% filter(Estimate=="I n,sco")
DataRes.n500.Iobs <- DataRes.n500 %>% filter(Estimate=="I n,obs")

mIsco500 <- round(c(mean(sqrt(500)*DataRes.n500.Isco$EstF11.true),mean(sqrt(500)*DataRes.n500.Isco$EstF22.true),
                    mean(sqrt(500)*DataRes.n500.Isco$EstF33.true),mean(sqrt(500)*DataRes.n500.Isco$EstF12.true),
                    mean(sqrt(500)*DataRes.n500.Isco$EstF13.true),mean(sqrt(500)*DataRes.n500.Isco$EstF23.true)),5)
mIobs500 <- round(c(mean(sqrt(500)*DataRes.n500.Iobs$EstF11.true),mean(sqrt(500)*DataRes.n500.Iobs$EstF22.true),
                    mean(sqrt(500)*DataRes.n500.Iobs$EstF33.true),mean(sqrt(500)*DataRes.n500.Iobs$EstF12.true),
                    mean(sqrt(500)*DataRes.n500.Iobs$EstF13.true),mean(sqrt(500)*DataRes.n500.Iobs$EstF23.true)),5)

sIsco500 <- round(c(sd(sqrt(500)*DataRes.n500.Isco$EstF11.true),sd(sqrt(500)*DataRes.n500.Isco$EstF22.true),
                    sd(sqrt(500)*DataRes.n500.Isco$EstF33.true),sd(sqrt(500)*DataRes.n500.Isco$EstF12.true),
                    sd(sqrt(500)*DataRes.n500.Isco$EstF13.true),sd(sqrt(500)*DataRes.n500.Isco$EstF23.true)),5)
sIobs500 <- round(c(sd(sqrt(500)*DataRes.n500.Iobs$EstF11.true),sd(sqrt(500)*DataRes.n500.Iobs$EstF22.true),
                    sd(sqrt(500)*DataRes.n500.Iobs$EstF33.true),sd(sqrt(500)*DataRes.n500.Iobs$EstF12.true),
                    sd(sqrt(500)*DataRes.n500.Iobs$EstF13.true),sd(sqrt(500)*DataRes.n500.Iobs$EstF23.true)),5)

dataBiasLMM <- cbind(component = c(1,2,3,4,5,6),
                     mIsco20 = mIsco20,
                     mIobs20 = mIobs20,
                     mIsco100 = mIsco100,
                     mIobs100 = mIobs100,
                     mIsco500 = mIsco500,
                     mIobs500 = mIobs500)

dataBiasLMM <- as.data.frame(dataBiasLMM)
dataBiasLMM$component <-factor(dataBiasLMM$component, levels=c(1,2,3,4,5,6),
                               labels=c("1"='(β,β)',
                                        "2"='(η2,η2)',
                                        "3"='(σ2,σ2)',
                                        "4"='(β,η2)',
                                        "5"='(β,σ2)',
                                        "6"='(η2,σ2)'))
dataBiasLMM <- flextable(dataBiasLMM, cwidth=1.2)
dataBiasLMM <- add_header_row(
  x = dataBiasLMM, values = c("", "n=20", "n=100", "n=500"),
  colwidths = c(1, 2, 2, 2))
dataBiasLMM <- set_header_labels(dataBiasLMM, mIsco20="Isco", mIobs20="Iobs", mIsco100="Isco", mIobs100="Iobs", mIsco500="Isco", mIobs500="Iobs", component="")
dataBiasLMM <- align(dataBiasLMM, part = "all", align = "center")
dataBiasLMM