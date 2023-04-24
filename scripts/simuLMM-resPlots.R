## Graphical representation of the empirical densities of the components of
## the FIM estimators computed in the linear mixed effects model while knowing 
## the true parameter values

source('functions/Fisher_LMM.R')
j     <- 12   
beta   <- 3
sigma2 <- 5 
eta2   <- 2

fisher <- Fisher_LMM(beta,sigma2,eta2,j) # computation of the exact FIM

load('Rfiles/simusLMM.Rdata')

DataRes.n <- DataRes %>% filter(n==500)

F22 <- ggplot(filter(DataRes.n), aes(sqrt(500)*(EstF22.true-fisher[2,2]), color=Estimate)) +
  geom_density(bw=0.3,alpha=0.6) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  xlim(-1.5,1.5) +
  ylim(0,2) +
  ggtitle(bquote('('~eta^2~','~eta^2~')')) +
  theme(legend.position = c(0,0.9), plot.title = element_text(size=20,face="bold"))

F33 <- ggplot(DataRes.n, aes(sqrt(500)*(EstF33.true-fisher[3,3]), color=Estimate)) +
  geom_density(bw=0.3,alpha=0.6) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  xlim(-1.5,1.5) +
  ylim(0,2) +
  ggtitle(bquote('('~sigma^2~','~ sigma^2~')')) +
  theme(legend.position = "none", plot.title = element_text(size=20,face="bold"))


F12 <- ggplot(DataRes.n, aes(sqrt(500)*(EstF12.true-fisher[1,2]), color=Estimate)) +
  geom_density(bw=0.3,alpha=0.6) +
  scale_fill_manual(values = c("#984EA3",'#E69F00')) +
  xlab("") +
  ylab("") +
  xlim(-1.5,1.5) +
  ylim(0,2) +
  ggtitle(bquote('('~beta~','~ eta^2~')')) +
  theme(legend.position = "none", plot.title = element_text(size=20,face="bold"))

plot_grid(F22, F33, F12, ncol = 3, nrow = 1)

