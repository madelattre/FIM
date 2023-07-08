dataSdLMM <- cbind(component = c(1,2,3,4,5,6),
                     sIsco20 = sIsco20,
                     sIobs20 = sIobs20,
                     sIsco100 = sIsco100,
                     sIobs100 = sIobs100,
                     sIsco500 = sIsco500,
                     sIobs500 = sIobs500)

dataSdLMM <- as.data.frame(dataSdLMM)
dataSdLMM$component <-factor(dataSdLMM$component, levels=c(1,2,3,4,5,6),
                               labels=c("1"='(β,β)',
                                        "2"='(η2,η2)',
                                        "3"='(σ2,σ2)',
                                        "4"='(β,η2)',
                                        "5"='(β,σ2)',
                                        "6"='(η2,σ2)'))
dataSdLMM <- flextable(dataSdLMM, cwidth=1.2)
dataSdLMM <- add_header_row(
  x = dataSdLMM, values = c("", "n=20", "n=100", "n=500"),
  colwidths = c(1, 2, 2, 2))
dataSdLMM <- set_header_labels(dataSdLMM, sIsco20="Isco", sIobs20="Iobs", sIsco100="Isco", sIobs100="Iobs", sIsco500="Isco", sIobs500="Iobs", component="")
dataSdLMM <- align(dataSdLMM, part = "all", align = "center")
dataSdLMM <- bold(dataSdLMM, j=c(2,4,6))
dataSdLMM <- width(dataSdLMM, width = 0.8)
dataSdLMM