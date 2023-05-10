load('Rfiles/ResNLMEexponentialBiasIsco.Rdata')
load('Rfiles/ResNLMEexponentialBiasIobs.Rdata')

nbiter <- 2000#dim(MbiasIsco)[3]
nbiterburnin <- 1000

DataResBias <- data.frame(
  BiasF11=c(abs(MbiasIsco[1,1,]),abs(MbiasIobs[1,1,])),
  BiasF22=c(abs(MbiasIsco[2,2,]),abs(MbiasIobs[2,2,])),
  BiasF33=c(abs(MbiasIsco[3,3,]),abs(MbiasIobs[3,3,])),
  BiasF44=c(abs(MbiasIsco[4,4,]),abs(MbiasIobs[4,4,])),
  BiasF55=c(abs(MbiasIsco[5,5,]),abs(MbiasIobs[5,5,])),
  BiasF66=c(abs(MbiasIsco[6,6,]),abs(MbiasIobs[6,6,])),
  Iter=c(seq(nbiterburnin+1,(nbiter+nbiterburnin),10),seq(nbiterburnin+1,nbiter+nbiterburnin,10)),
  Estimate=c(rep('Isco',nbiter),rep('Iobs',nbiter)))

BiasF11 <- ggplot(DataResBias, aes(y=BiasF11, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~ka~','~ ka~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

BiasF22 <- ggplot(DataResBias, aes(y=BiasF22, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~V~','~ V~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

BiasF33 <- ggplot(DataResBias, aes(y=BiasF33, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~Cl~','~ Cl~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

BiasF44 <- ggplot(DataResBias, aes(y=BiasF44, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~omega[ka]^2~','~ omega[ka]^2~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

BiasF55 <- ggplot(DataResBias, aes(y=BiasF55, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~omega[V]^2~','~ omega[V]^2~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

BiasF66 <- ggplot(DataResBias, aes(y=BiasF66, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~omega[Cl]^2~','~ omega[Cl]^2~')')) +
  theme(plot.title = element_text(size=20,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))


plot_grid(BiasF11, BiasF22, BiasF33, BiasF44, BiasF55, BiasF66, ncol = 3, nrow = 2)

