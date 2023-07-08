dataSdMixt <- cbind(component = c(1,2,3,4,5,6),
                      sIsco20 = c(diag(sdIscoMixt20)[2:5],sdIscoMixt20[2,3],sdIscoMixt20[3,5]),
                      sIobs20 = c(diag(sdIobsMixt20)[2:5],sdIobsMixt20[2,3],sdIobsMixt20[3,5]),
                      sIsco100 =c(diag(sdIscoMixt100)[2:5],sdIscoMixt100[2,3],sdIscoMixt100[3,5]),
                      sIobs100 =c(diag(sdIobsMixt100)[2:5],sdIobsMixt100[2,3],sdIobsMixt100[3,5]),
                      sIsco500 =c(diag(sdIscoMixt500)[2:5],sdIscoMixt500[2,3],sdIscoMixt500[3,5]),
                      sIobs500 =c(diag(sdIobsMixt500)[2:5],sdIobsMixt500[2,3],sdIobsMixt500[3,5]))


dataSdMixt <- as.data.frame(dataSdMixt)
dataSdMixt$component <-factor(dataSdMixt$component, levels=c(1,2,3,4,5,6),
                                labels=c("1"='(λ2,λ2)',
                                         "2"='(λ3,λ3)',
                                         "3"='(α1,α1)',
                                         "4"='(α2,α2)',
                                         "5"='(λ2,λ3)',
                                         "6"='(λ3,α2)'))
dataSdMixt <- flextable(dataSdMixt, cwidth=1.2)
dataSdMixt <- add_header_row(
  x = dataSdMixt, values = c("", "n=20", "n=100", "n=500"),
  colwidths = c(1, 2, 2, 2))
dataSdMixt <- set_header_labels(dataSdMixt, sIsco20="Isco", sIobs20="Iobs", 
                                  sIsco100="Isco", sIobs100="Iobs", sIsco500="Isco", 
                                  sIobs500="Iobs", component="")
dataSdMixt <- align(dataSdMixt, part = "all", align = "center")
dataSdMixt <- bold(dataSdMixt, j=c(2,4,6))
dataSdMixt <- width(dataSdMixt, width = 0.8)
dataSdMixt