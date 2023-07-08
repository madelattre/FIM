## Nonlinear mixed-effects model not belonging to the curved exponential family
## Convergence graphs

load("Rfiles/saem-non-exp-conv-plot.Rdata")

nbiterem <- 3000

DataResEst <- data.frame(
  ka       = res$thetaest[1,],
  V        = res$thetaest[2,],  
  Cl       = res$thetaest[3,],
  o2ka     = res$thetaest[4,],
  Iscoka   = res$isco[1,1,],
  Iscov    = res$isco[2,2,],
  Iscocl   = res$isco[3,3,],
  Iscoo2ka = res$isco[4,4,],
  Iter = c(1,seq(10,nbiterem,10)))

kaconv <- ggplot(DataResEst, aes(y=ka, x=Iter)) +
  geom_line(size=0.5) +
  xlab("") +
  ylab("") +
  ggtitle(bquote(ka)) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=45))

vconv <- ggplot(DataResEst, aes(y=V, x=Iter)) +
  geom_line(size=0.5) +
  xlab("") +
  ylab("") +
  ggtitle(bquote(V)) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=45))

clconv <- ggplot(DataResEst, aes(y=Cl, x=Iter)) +
  geom_line(size=0.5) +
  xlab("") +
  ylab("") +
  ggtitle(bquote(Cl)) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=45))

o2kaconv <- ggplot(DataResEst, aes(y=o2ka, x=Iter)) +
  geom_line(size=0.5) +
  xlab("") +
  ylab("") +
  ggtitle(bquote(omega[ka]^2)) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=45))

iscokaconv <- ggplot(DataResEst, aes(y=Iscoka, x=Iter)) +
  geom_line(size=0.5) +
  xlab("") +
  ylab("") +
  ggtitle(bquote(I[n-sco](ka,ka))) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=45))

iscovconv <- ggplot(DataResEst, aes(y=Iscov, x=Iter)) +
  geom_line(size=0.5) +
  xlab("") +
  ylab("") +
  #ylim(0.1,0.25) +
  ggtitle(bquote(I[n-sco](V,V))) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=45))

iscoclconv <- ggplot(DataResEst, aes(y=Iscocl, x=Iter)) +
  geom_line(size=0.5) +
  xlab("") +
  ylab("") +
  ggtitle(bquote(I[n-sco](Cl,Cl))) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=45))

iscoo2kaconv <- ggplot(DataResEst, aes(y=Iscoo2ka, x=Iter)) +
  geom_line(size=0.5) +
  xlab("") +
  ylab("") +
  ggtitle(bquote(I[n-sco](omega[ka]^2,omega[ka]^2))) +
  theme(plot.title = element_text(size=10,face="bold"),legend.position='none',
        axis.text.x = element_text(angle=45))

plot_grid(kaconv, vconv, clconv, o2kaconv, iscokaconv, iscovconv, iscoclconv, 
          iscoo2kaconv, ncol = 4, nrow = 2)

