library(phytools)
library(fftw)
set.seed(7996)

k<-2
n<-50
tree<-pbtree(n=n,scale=1)
Q<-rand.Q(k,5,5)
sig2<-rgamma(k,0.2,0.2)
sim<-sim.sce(tree,Q,sig2)
plot(sim,ftype='off')

sce<-make.sce(tree,sim[['disc']],sim[['cont']],hidden=TRUE,tip.sig=0.1)
fit<-find.mle(sce,control=list(trace=TRUE))
recon<-anc.recon(fit)
phenogram(sim$simmap,recon$anc.recon$means[,k+1],ftype='off')
nodelabels(pie=recon$disc.recon[-seq_len(ntips),])
plot(sim$simmap)
tiplabels(pie=recon$disc.recon[seq_len(ntips),],cex=0.5,piecol=c('red','black')) #so cool!
