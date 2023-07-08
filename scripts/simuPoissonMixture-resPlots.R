## Graphical representation of the empirical distribution of the FIM estimates
## in the Poisson mixture model

n <- 500
filename <- paste('Rfiles/simusMixt_n',n,'.Rdata',sep="")
load(filename)
Isco <- ResSim$Isco
Iobs <- ResSim$Iobs
load('Rfiles/PoissonMixtureTrueFIM.Rdata')
nbsim <- length(Isco[1,1,])

DataRes <- data.frame(EstF11=c(sqrt(n)*(Isco[1,1,]-trueFIM[1,1]),
                               sqrt(n)*(Iobs[1,1,]-trueFIM[1,1])),
                      EstF22=c(sqrt(n)*(Isco[2,2,]-trueFIM[2,2]),
                               sqrt(n)*(Iobs[2,2,]-trueFIM[2,2])),
                      EstF33=c(sqrt(n)*(Isco[3,3,]-trueFIM[3,3]),
                               sqrt(n)*(Iobs[3,3,]-trueFIM[3,3])),
                      EstF44=c(sqrt(n)*(Isco[4,4,]-trueFIM[4,4]),
                               sqrt(n)*(Iobs[4,4,]-trueFIM[4,4])),
                      EstF55=c(sqrt(n)*(Isco[5,5,]-trueFIM[5,5]),
                               sqrt(n)*(Iobs[5,5,]-trueFIM[5,5])),
                      EstF12=c(sqrt(n)*(Isco[1,2,]-trueFIM[1,2]),
                               sqrt(n)*(Iobs[1,2,]-trueFIM[1,2])),
                      EstF13=c(sqrt(n)*(Isco[1,3,]-trueFIM[1,3]),
                               sqrt(n)*(Iobs[1,3,]-trueFIM[1,3])),
                      EstF23=c(sqrt(n)*(Isco[2,3,]-trueFIM[2,3]),
                               sqrt(n)*(Iobs[2,3,]-trueFIM[2,3])),
                      EstF35=c(sqrt(n)*(Isco[3,5,]-trueFIM[3,5]),
                               sqrt(n)*(Iobs[3,5,]-trueFIM[3,5])),
                      EstF34=c(sqrt(n)*(Isco[3,4,]-trueFIM[3,4]),
                               sqrt(n)*(Iobs[3,4,]-trueFIM[3,4])),
                      EstF25=c(sqrt(n)*(Isco[2,5,]-trueFIM[2,5]),
                               sqrt(n)*(Iobs[2,5,]-trueFIM[2,5])),
                      EstF24=c(sqrt(n)*(Isco[2,4,]-trueFIM[2,4]),
                               sqrt(n)*(Iobs[2,4,]-trueFIM[2,4])),
                      EstF15=c(sqrt(n)*(Isco[1,5,]-trueFIM[1,5]),
                               sqrt(n)*(Iobs[1,5,]-trueFIM[1,5])),
                      EstF14=c(sqrt(n)*(Isco[1,4,]-trueFIM[1,4]),
                               sqrt(n)*(Iobs[1,4,]-trueFIM[1,4])),
                      Estimate=c(rep('I N,sco',nbsim),rep('I N,obs',nbsim)))


F11 <- ggplot(DataRes, aes(EstF11, color=Estimate)) +
  geom_density(bw=0.2,linewidth=0.6,
               key_glyph = draw_key_path) +
  scale_color_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  xlim(-1,1) +
  ylim(0,2.3) +
  ggtitle(bquote('('~lambda[1]~','~ lambda[1]~')')) +
  theme(legend.position = c(0,1),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        plot.title = element_text(size=10,face="bold"))

F22 <- ggplot(DataRes, aes(EstF22, color=Estimate)) +
  geom_density(bw=0.2,linewidth=0.6,
               key_glyph = draw_key_path) +
  scale_color_manual(values = c("#984EA3",'#E69F00'))  +
  xlab("") +
  ylab("") +
  xlim(-1,1) +
  ylim(0,2.3) +
  ggtitle(bquote('('~lambda[2]~','~ lambda[2]~')')) +
  theme(legend.position = "none", 
        plot.title = element_text(size=10,face="bold"))

F25 <- ggplot(DataRes, aes(EstF25, color=Estimate)) +
  geom_density(bw=0.2,linewidth=0.6,
               key_glyph = draw_key_path) +
  scale_color_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  xlim(-1.5,1.5) +
  ylim(0,2.3) +
  ggtitle(bquote('('~lambda[2]~','~ alpha[2]~')')) +
  theme(legend.position = "none", 
        plot.title = element_text(size=10,face="bold"))


plot_grid(F11, F22, F25, ncol = 3, nrow = 1)

