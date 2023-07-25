rm(list=ls())
library(phytools)
library(sce)
tree<-pbtree(n=50,scale=1)
test<-sim.sce(tree,Q=2*matrix(c(-1,1,1,-1),2,2),sig2=c(1,20))
plot(test,ftype='off')

library(expm)
fxn<-make.sce2(tree,test$disc,test$cont,nx=128)
test1<-optim(rnorm(4),fxn$lik.fun,method='BFGS',control=list(trace=TRUE)) #kinda slow, but overall not too bad

library(fftw)
fxn2<-make.sce(tree,test$disc,test$cont,nx=1024,nt=500)
test2<-optim(rnorm(4),fxn2$lik.fun,control=list(trace=TRUE)) #it's actually really close to the true value...wow!

#why is expokit so slow? Must have to use it in C directly?
library(rexpokit)
Q<-matrix(rexp(16),4,4)
diag(Q)<-0
diag(Q)<- -rowSums(Q)
get.Q<-get("get.Q",env=environment(fxn$lik.fun))
tmp<-mat2coo(get.Q(log(c(1,2,2,20))))
expokit_dgexpv_Qmat(tmp,t=1,transform_to_coo_TF=FALSE,coo_n=256,anorm=1)

#alright...yeah--at least you conclusively demonstrated to yourself that your approximation method is legit
#virtually the same parameters, at higher resolution, in a fraction of the time
fxn2<-make.sce(tree,test$disc,test$cont)
test<-find.mle(fxn2,control=list(trace=TRUE))
testy<-anc.recon(test)
