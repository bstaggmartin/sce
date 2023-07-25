library(phytools)
library(fftw)
library(sce)
library(ucminf)
set.seed(123)
tree<-pbtree(n=100,scale=1)
par<-c(2,10,0)
dat<-sim.levy(tree,model=2,rate=par[1],jump=par[2],sig2=par[3])
phenogram(tree,dat[,1],ftype='off')
levy<-make.levy(tree,dat[,1],model=2,nx=1024,tip.sig=0.1)
levy[['lik.fun']](par)
# tmp<-optim(rexp(3),levy[['lik.fun']],control=list(trace=TRUE),method='BFGS',hessian=TRUE)
tmp<-ucminf(rexp(3),levy[['lik.fun']],control=list(trace=TRUE),hessian=TRUE)
(tmp$par)
(cbind(tmp$par[-3]-2*sqrt(diag(solve(tmp$hessian[-3,-3]))),
          tmp$par[-3]+2*sqrt(diag(solve(tmp$hessian[-3,-3]))))) #hessian seems to indicate high uncertainty, so the bias is probably okay
#bias has something to do with estimating variance params, probably
#definitely need to check DFTs against simulations to ensure proper scaling, but it actually looks pretty good
#the likelihoods seem to be on very different scales across levy models, unsure what that's about
test<-levy[['recon.fun']](tmp$par)
xpts<-attr(test,'xpts')
dx<-diff(xpts[1:2])
anc.means<-dx*.rowSums(sweep(test,2,xpts,'*'),dim(test)[1],dim(test)[2])
nodes<-get('nodes',envir=environment(levy[['recon.fun']]),inherits=FALSE)[,2]
ord<-order(nodes)
anc.means<-anc.means[ord]
names(anc.means)<-c(tree[['tip.label']],seq_len(tree[['Nnode']])+Ntip(tree))
rbind(anc.means[tree$tip.label],dat[,1]) #looks right
phenogram(tree,anc.means,ftype='off')
tips<-get('tips',envir=environment(levy[['recon.fun']]),inherits=FALSE)
matplot(t(test[-tips,]),type='l',x=xpts) #some bimodality for compound poisson! cool!

test<-levy[['map.fun']](tmp$par,nt=100,stochastic=TRUE,nsim=20)
des<-get('des',envir=environment(levy[['map.fun']]))
sample.des.path<-function(des){
  out<-NULL
  tmp.des<-des[[length(des)]]
  while(length(tmp.des)){
    cur.edge<-sample(tmp.des,1)
    out<-c(out,cur.edge)
    tmp.des<-des[[cur.edge]]
  }
  out
}
viz.pdf.path<-function(max.prob=0.1){
  path<-sample.des.path(des)
  pdf.path<-test[path]
  tt<-lapply(pdf.path,function(ii) attr(ii,'tpts'))
  lens<-lengths(tt)
  nedge<-length(lens)
  tmp.seq<-seq_len(nedge-1)
  tt<-c(unlist(lapply(tmp.seq,function(ii) tt[[ii]][-lens[ii]])),tt[[nedge]])
  xx<-attr(test,'xpts')
  pdf.path<-do.call(rbind,c(lapply(tmp.seq,function(ii) pdf.path[[ii]][-lens[ii],,drop=FALSE]),
                            list(pdf.path[[nedge]])))
  breaks<-c(-1,seq(0,max.prob/diff(xpts[1:2]),length.out=100),1/diff(xpts[1:2]))
  col<-c('black',colorRampPalette(c('black','white'))(99),'white')
  image(pdf.path,y=xx,x=tt,breaks=breaks,col=col,
        xlab='time',
        ylab='trait value')
  legend('topright',fill=c('black','white'),legend=c(' = 0',paste(' >',max.prob)),
         border='white',text.col='white',title.col='white',bty='n',
         title='posterior probability')
}
viz.pdf.path(max.prob=0.1)

xpts<-attr(test,'xpts')
dx<-diff(xpts[1:2])
means<-lapply(test,function(ii) dx*rowSums(sweep(ii,2,xpts,'*')))
tt<-lapply(test,function(ii) attr(ii,'tpts'))

plot(0,col='white',ylim=range(unlist(means)),xlim=range(unlist(tt)))
for(i in seq_along(means)){
  lines(as.vector(means[[i]])~tt[[i]],col='black')
}

cdfs<-lapply(test,function(ii) t(apply(ii,1,cumsum)*dx))
matplot(t(cdfs[[50]]),type='l')
get.quant<-function(prob){
  lapply(cdfs,function(ii) t(apply(ii,1,function(jj) xpts[which.min(abs(prob-jj))])))
}
means<-get.quant(0.5)

simmaps<-test
tt<-lapply(simmaps,function(ii) attr(ii,'tpts'))
plot(0,col='white',ylim=range(unlist(simmaps)),xlim=range(unlist(tt)),
     xlab='time',ylab='trait value')
for(i in seq_along(simmaps)){
  matplot(simmaps[[i]][,seq_len(20),drop=FALSE],x=tt[[i]],col=rgb(0,0,0,0.1),add=TRUE,lty=1,type='l')
}
matplot(sim.levy.path(rate=tmp$par[1],jump=tmp$par[2],sig2=tmp$par[3],model=2,nsim=25),type='l',lty=1,col='black')
