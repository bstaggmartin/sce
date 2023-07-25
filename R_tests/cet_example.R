library(evorates)
library(sce)
library(phytools)
library(fftw)

data("cet_fit")
tree<-cet_fit$call$tree
dat<-setNames(cet_fit$call$trait.data,
              rownames(cet_fit$call$trait.data))

set.seed(123456)
Q<-rand.Q(2,5,100)
rownames(Q)<-colnames(Q)<-c("0","1")
map<-sim.history(tree,Q)
phenogram(tree,dat,ftype="off",
          ylab="log body size (m)")
tiplabels(pch=16,col=setNames(c("black","red"),c("0","1"))[map$states],
          cex=2,
          adj=0)
tiplabels(pie=rep(0.5,length(tree$tip.label)),piecol=c("green","blue"),
          cex=0.8,
          adj=1.2)

tmp<-make.sce(tree,map$states,dat,tip.sig=0.05,nx=256,
              Q.template=matrix(c(1,2,3,1),2,2))
fit0<-find.mle(tmp,inits.mean=0.05,control=list(trace=TRUE))
tmp<-make.sce(tree,map$states,dat,tip.sig=0.05)
fit1<-find.mle(tmp,inits.mean=0.05,control=list(trace=TRUE))
AIC<-function(...,small.sample=TRUE){
  nms<-strsplit(gsub("\\)$","",gsub("^list\\(","",deparse(substitute(list(...))))),", ")[[1]]
  tmp<-list(...)
  lnL<-sapply(tmp,'[[',"lnL")
  k<-sapply(tmp,function(ii) length(ii$par))
  n<-sapply(tmp,function(ii) length(get("tips",environment(ii$sce$lik.fun))))
  sort(setNames(2*(k-lnL+k*(k+1)/(n-k-1)),nms))
}
AIC(fit1,fit0) #fit1 "significantly" better
states<-map$states
states<-setNames(c("0A&0B","1A&1B")[(states=="1")+1],names(states))
Q.template<-matrix(c(1,2,3,0,
                     4,5,0,3,
                     6,0,7,2,
                     0,6,4,8),
                   4,4)
rownames(Q.template)<-colnames(Q.template)<-c("0A","1A","0B","1B")
tmp<-make.sce(tree,states,dat,tip.sig=0.05,nx=256,
              Q.template=Q.template)
fit2<-find.mle(tmp,inits.mean=0.05,control=list(trace=TRUE))
Q.template<-matrix(c(1,2,3,0,
                     4,5,0,3,
                     6,0,1,2,
                     0,6,4,5),
                   4,4)
rownames(Q.template)<-colnames(Q.template)<-c("0A","1A","0B","1B")
tmp<-make.sce(tree,states,dat,tip.sig=0.05,nx=256,
              Q.template=Q.template)
fit3<-find.mle(tmp,inits.mean=0.05,control=list(trace=TRUE))
Q.template<-matrix(c(1,2,3,0,
                     4,1,0,3,
                     5,0,6,2,
                     0,5,4,6),
                   4,4)
rownames(Q.template)<-colnames(Q.template)<-c("0A","1A","0B","1B")
tmp<-make.sce(tree,states,dat,tip.sig=0.05,nx=256,
              Q.template=Q.template)
fit4<-find.mle(tmp,inits.mean=0.05,control=list(trace=TRUE))
AIC(fit4,fit3,fit2)

recon<-anc.recon(fit4)
disc.recon<-cbind(rowSums(recon$disc.recon[,c("0A","1A")]),
                  rowSums(recon$disc.recon[,c("0B","1B")]))
colnames(disc.recon)<-c("A","B")
cont.recon<-recon$anc.recon$means[,"marginal"]
phenogram(tree,cont.recon,ftype="off",ylab="log body size (m)")
tmp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cont.vars<-recon$anc.recon$vars[,"marginal"]
segments(x0=tmp$xx,y0=tmp$yy-2*sqrt(cont.vars),y1=tmp$yy+2*sqrt(cont.vars),lty=2)
nodelabels(pie=disc.recon[-seq_along(tree$tip.label),],piecol=c("blue","green"))
tiplabels(pie=disc.recon[seq_along(tree$tip.label),],piecol=c("blue","green"))
plot(ladderize(tree),show.tip.label=FALSE)
nodelabels(pie=disc.recon[-seq_along(tree$tip.label),],piecol=c("blue","green"))
tiplabels(pie=disc.recon[seq_along(tree$tip.label),],piecol=c("blue","green"),cex=0.5)

#further analysis REALLY seems to support that there are only 2 rates on this tree! Crazy...
#maybe try without tip.sig?
#yeah--without tip.sig, support for three rates, with delphinids showing equivocal support for somewhat higher rates...
Q.template<-matrix(c(1,2,0,
                     3,4,2,
                     0,3,5),
                   3,3)
rownames(Q.template)<-colnames(Q.template)<-c("A","B","C")
tmp.states<-states
tmp.states[]<-c("A&B&C")
tmp<-make.sce(tree,tmp.states,dat,tip.sig=0,nx=256,
              Q.template=Q.template)
fit5<-find.mle(tmp,inits.mean=0.05,control=list(trace=TRUE))
recon<-anc.recon(fit5)
cont.recon<-recon$anc.recon$means[,"marginal"]
phenogram(tree,cont.recon,ftype="off",ylab="log body size (m)")
tmp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cont.vars<-recon$anc.recon$vars[,"marginal"]
segments(x0=tmp$xx,y0=tmp$yy-2*sqrt(cont.vars),y1=tmp$yy+2*sqrt(cont.vars),lty=2)
plot(ladderize(tree),show.tip.label=FALSE)
nodelabels(pie=recon$disc.recon[-seq_along(tree$tip.label),],piecol=c("red","blue","green"))
tiplabels(pie=recon$disc.recon[seq_along(tree$tip.label),],piecol=c("red","blue","green"))
