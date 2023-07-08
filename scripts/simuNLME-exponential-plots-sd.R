load('Rfiles/ResNLMEexponentialSdIsco.Rdata')
load('Rfiles/ResNLMEexponentialSdIobs.Rdata')

DataResRsd <- data.frame(
  RsdF11=c(abs(MsdIsco[1,1,-seq(1,nbiterburnin)]),abs(MsdIobs[1,1,-seq(1,nbiterburnin)])),
  RsdF22=c(abs(MsdIsco[2,2,-seq(1,nbiterburnin)]),abs(MsdIobs[2,2,-seq(1,nbiterburnin)])),
  RsdF33=c(abs(MsdIsco[3,3,-seq(1,nbiterburnin)]),abs(MsdIobs[3,3,-seq(1,nbiterburnin)])),
  RsdF44=c(abs(MsdIsco[4,4,-seq(1,nbiterburnin)]),abs(MsdIobs[4,4,-seq(1,nbiterburnin)])),
  RsdF55=c(abs(MsdIsco[5,5,-seq(1,nbiterburnin)]),abs(MsdIobs[5,5,-seq(1,nbiterburnin)])),
  RsdF66=c(abs(MsdIsco[6,6,-seq(1,nbiterburnin)]),abs(MsdIobs[6,6,-seq(1,nbiterburnin)])),
  Iter=c(seq(nbiterburnin+1,nbiter,1),seq(nbiterburnin+1,nbiter,1)),
  Estimate=c(rep('Isco',nbiter-nbiterburnin),rep('Iobs',nbiter-nbiterburnin)))


RsdF11 <- ggplot(DataResRsd, aes(y=RsdF11, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~ka~','~ ka~')')) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

RsdF22 <- ggplot(DataResRsd, aes(y=RsdF22, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~V~','~ V~')')) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

RsdF33 <- ggplot(DataResRsd, aes(y=RsdF33, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~Cl~','~ Cl~')')) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

RsdF44 <- ggplot(DataResRsd, aes(y=RsdF44, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~omega[ka]^2~','~ omega[ka]^2~')')) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

RsdF55 <- ggplot(DataResRsd, aes(y=RsdF55, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~omega[V]^2~','~ omega[V]^2~')')) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

RsdF66 <- ggplot(DataResRsd, aes(y=RsdF66, x=Iter, color=Estimate)) +
  geom_line(size=1) +
  xlab("") +
  ylab("") +
  ggtitle(bquote('('~omega[Cl]^2~','~ omega[Cl]^2~')')) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=90))

plot_grid(RsdF11, RsdF22, RsdF33, RsdF44, RsdF55, RsdF66, ncol = 3, nrow = 2)
