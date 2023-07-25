library(phytools)
library(fftw)
library(sce)
set.seed(7996)
# set.seed(123)

k<-2
n<-50
tree<-pbtree(n=n,scale=1)
Q<-rand.Q(k,5,5)
sig2<-rgamma(k,0.5,0.5)
sim<-sim.sce(tree,Q,sig2)
plot(sim,ftype='off')
# disc<-setNames(c('1A&1B','2A&2B')[as.numeric(sim[['disc']])],names(sim[['disc']]))
# Q.template<-matrix(c(1,2,3,0,4,5,0,3,6,0,7,2,0,6,4,8),2*k,2*k,
#                    dimnames=rep(list(c('1A','2A','1B','2B')),2))
disc2<-setNames(rep(c('1&2'),n),names(sim[['disc']]))
sce<-make.sce(tree,disc2,sim[['cont']],nx=1024,tip.sig=0.05)
tmp<-Q
diag(tmp)<-sig2
par<-as.vector(log(tmp))
sce$lik.fun(par)
fit<-find.mle(sce,hessian=FALSE,control=list(trace=TRUE),method="BFGS")
fit<-find.mle(sce,hessian=FALSE,control=list(trace=TRUE),method="Nelder-Mead")

debug(sce$lik.fun)
fit<-find.mle(sce,hessian=FALSE,control=list(trace=TRUE),method="BFGS")
#cranking down nx seems to preserve relative likelihoods, so you can still get a decent ballpark estimate
#but the approximation to the true likelihood is pretty bad!
#still, suggests a workflow for high-dimensional problems where you first find the "support set" using a low-dimensional approximation, then
##find the "true likelihood" calculating at high dimensions
#also, likelihoods are much more numerically stable with non-0 tip.sig (makes sense given rounding...)
#may be better to automatically use a tip.sig proportional to dx when it's set to 0? Unsure...
#a little over a minute to find MLE with nx=1024 (50 tips, 2 discrete states, using BFGS)
fit<-list('sce'=sce,'par'=par)
recon<-anc.recon(fit)
k<-length(attr(sce,'key'))
phenogram(sim$simmap,recon$anc.recon$means[,k+1],ftype='off')
segments(x0=node.depth.edgelength(tree),
         y0=recon$anc.recon$means[,k+1]-2*sqrt(recon$anc.recon$vars[,k+1]),
         y1=recon$anc.recon$means[,k+1]+2*sqrt(recon$anc.recon$vars[,k+1]))
nodelabels(pie=recon$disc.recon[-seq_len(n),],piecol=seq_len(k),cex=0.4)
tiplabels(pie=recon$disc.recon[seq_len(n),],cex=0.5,piecol=seq_len(k),cex=0.2)

pdf('~/../Desktop/example_phy.pdf',height=10,width=10)
plot(sim$simmap,ftype='off')
nodelabels(pie=recon$disc.recon[-seq_len(n),],piecol=seq_len(k),cex=0.5)
tiplabels(pie=recon$disc.recon[seq_len(n),],piecol=seq_len(k),cex=0.5) #so cool!
dev.off()
pdf('~/../Desktop/example_phen.pdf',height=10,width=10)
phenogram(tree,sim$cont,ftype='off')
dev.off()

tmp<-recon$ancestral.pdfs$y[,1,]
matplot(tmp[,-seq_len(n)],x=recon$ancestral.pdfs$x,type='l',add=TRUE)
matplot(log(tmp)[,-seq_len(n)],x=recon$ancestral.pdfs$x,type='l')
matplot(tmp,x=recon$ancestral.pdfs$x,type='l')


####NEW IDEA THAT WILL TAKE SOME REWORKING OF THINGS####
#initialize sce's with "low res" and "hi res" versions of likelihood function
#find.mle() will then use the "low res" version with random restarts to get a good initial guess
#the "high res" version will be called to get the final guess
#initial tests seem to indicate that this would be a decent strategy!

#might not have fully figured out how to normalize for hidden rates yet...but it could also be numerical error?
#yeah...I think it's some kind of numerical issue--probably rerunning some more times it would find a more comparable solution
#oh, maybe not...it does have to explain the evolution of the discrete trait data too now, so that's probably why
#you could compare it to a no hidden states setup, but not a hidden states only setup (though you could constrain sig2 parameters to get at that)
