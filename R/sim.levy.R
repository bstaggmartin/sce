.rinvgauss<-function(n,mean,shape){
  nu<-rnorm(n)
  nu.sq<-nu^2
  mean.sq<-mean^2
  meannu.sq<-mean.sq*nu.sq
  shape2<-shape*2
  x<-mean+meannu.sq/shape2-mean/shape2*sqrt(4*mean*shape*nu.sq+meannu.sq*nu.sq)
  z<-runif(n)
  inds<-z>mean/(mean+x)
  x[inds]<-rep(mean.sq,length.out=length(x))[inds]/x[inds]
  x
}

.get.seed<-function(model,
                    rate,jump,sig2,root,
                    dt,nt,nsim){
  foo<-function(x) rep(x,each=nt)
  rate<-foo(rate);jump<-foo(jump);sig2<-foo(sig2)
  scalars<-sqrt(switch(model,
                       jump*rpois(nt*nsim,dt*rate)+dt*sig2,
                       .rinvgauss(nt*nsim,dt*rate*sqrt(jump),(dt*rate)^2)+dt*sig2,
                       jump*rgamma(nt*nsim,dt*rate,1/rate)+dt*sig2))
  seed<-matrix(scalars*rnorm(nt*nsim),nt,nsim)
  seed
}

#rate roughly controls jump frequency
#jump roughly controls jump magnitude/shape
#this breaks down a bit in NIG, where rate corresponds to the scale of small jumps and jump to scale of large jumps
#the frequency of jumps under NIG seems rather constant...
#rate = lambda for JN, delta for NIG, and 1/nu for VG
#jump = tau2 for JN, 1/alpha2 for NIG, and tau2 for VG
#' @export
sim.levy<-function(tree,rate=1,jump=1,sig2=0,root=0,
                   model=c('JN','NIG','VG'),internal=FALSE,nsim=1){
  list2env(.get.tree.topo.info(tree),envir=environment()) #gets more info than necessary...only need ancestor indices...
  choices<-c('JN','NIG','VG')
  if(is.character(model)){
    model<-pmatch(model[1],choices)
  }
  seed<-.get.seed(model,
                  rate,jump,sig2,root,
                  elen,nedge-1,nsim)
  out<-matrix(nrow=nedge,ncol=nsim)
  out[nedge,]<-root
  for(i in clade.seq){
    out[i,]<-out[anc[i],,drop=FALSE]+seed[i,,drop=FALSE]
  }
  inds<-nodes[,2,drop=FALSE]
  ord<-order(inds)
  inds[tips]<-tips.ord
  rownames(out)<-inds
  out<-out[ord,,drop=FALSE]
  if(internal){
    out
  }else{
    out[seq_len(length(tips)),,drop=FALSE]
  }
}

#' @export
sim.levy.path<-function(t=1,res=100,rate=1,jump=1,sig2=0,x0=0,
                        model=c('JN','NIG','VG'),nsim=1){
  tt<-seq(0,t,length.out=res)
  dt<-diff(tt[1:2])
  choices<-c('JN','NIG','VG')
  if(is.character(model)){
    model<-pmatch(model[1],choices)
  }
  seed<-.get.seed(model,
                  rate,jump,sig2,root,
                  dt,res-1,nsim)
  out<-rbind(x0,seed)
  apply(out,2,cumsum)
}
