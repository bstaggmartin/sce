#' @export
rand.Q<-function(k,shape,rate){
  qq<-rgamma(k^2-k,shape,rate)
  Q<-matrix(0,k,k)
  Q.diag<-seq.int(1,k^2,k+1)
  Q.inds<-seq_len(k^2)[-Q.diag]
  Q[Q.inds]<-qq
  Q[Q.diag]<- -.rowSums(Q,k,k)
  Q
}

#' @export
sim.sce<-function(tree,Q,sig2,disc.anc=NULL,cont.anc=0){
  simmap<-phytools::sim.history(tree,Q,anc=disc.anc)
  disc<-as.factor(simmap$states)
  if(is.null(names(sig2))){
    names(sig2)<-levels(disc)
  }
  cont<-phytools::sim.rates(simmap,sig2,anc=cont.anc)
  out<-list('disc'=disc,'cont'=cont,'simmap'=simmap)
  class(out)<-'sce_sim'
  out
}

#' @export
plot.sce_sim<-function(sim,...){
  phenogram(sim[['simmap']],sim[['cont']],...)
}