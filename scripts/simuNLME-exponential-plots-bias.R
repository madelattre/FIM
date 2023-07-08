load('Rfiles/ResNLMEexponentialBiasIsco.Rdata')
load('Rfiles/ResNLMEexponentialBiasIobs.Rdata')

nbiter <- dim(MbiasIsco)[3]
nbiterburnin <- 1000

DataResBias <- data.frame(
  BiasF11=c(abs(MbiasIsco[1,1,-seq(1,nbiterburnin)]),abs(MbiasIobs[1,1,-seq(1,nbiterburnin)])),
  BiasF22=c(abs(MbiasIsco[2,2,-seq(1,nbiterburnin)]),abs(MbiasIobs[2,2,-seq(1,nbiterburnin)])),
  BiasF33=c(abs(MbiasIsco[3,3,-seq(1,nbiterburnin)]),abs(MbiasIobs[3,3,-seq(1,nbiterburnin)])),
  BiasF44=c(abs(MbiasIsco[4,4,-seq(1,nbiterburnin)]),abs(MbiasIobs[4,4,-seq(1,nbiterburnin)])),
  BiasF55=c(abs(MbiasIsco[5,5,-seq(1,nbiterburnin)]),abs(MbiasIobs[5,5,-seq(1,nbiterburnin)])),
  BiasF66=c(abs(MbiasIsco[6,6,-seq(1,nbiterburnin)]),abs(MbiasIobs[6,6,-seq(1,nbiterburnin)])),
  Iter=c(seq(nbiterburnin+1,nbiter,1),seq(nbiterburnin+1,nbiter,1)),
  Estimate=c(rep('Isco',nbiter-nbiterburnin),rep('Iobs',nbiter-nbiterburnin)))

BiasF11 <- ggplot(DataResBias, aes(y=BiasF11, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~ka~','~ ka~')')) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

BiasF22 <- ggplot(DataResBias, aes(y=BiasF22, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~V~','~ V~')')) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

BiasF33 <- ggplot(DataResBias, aes(y=BiasF33, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~Cl~','~ Cl~')')) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

BiasF44 <- ggplot(DataResBias, aes(y=BiasF44, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~omega[ka]^2~','~ omega[ka]^2~')')) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

BiasF55 <- ggplot(DataResBias, aes(y=BiasF55, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~omega[V]^2~','~ omega[V]^2~')')) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

BiasF66 <- ggplot(DataResBias, aes(y=BiasF66, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~omega[Cl]^2~','~ omega[Cl]^2~')')) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))


plot_grid(BiasF11, BiasF22, BiasF33, BiasF44, BiasF55, BiasF66, ncol = 3, nrow = 2)

