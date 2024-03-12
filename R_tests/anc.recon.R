#' @export
anc.recon<-function(fit,par=NULL){
  key<-attr(fit[['sce']],'key')
  k<-length(key)
  recon.fun<-fit[['sce']][['recon.fun']]
  if(is.null(par)){
    par<-c(fit[['par']])
  }
  PP<-recon.fun(par)
  xpts<-attr(PP,'xpts')
  xpts.sq<-xpts^2
  dx<-diff(xpts[1:2])

  tree<-get('tree',envir=environment(recon.fun),inherits=FALSE)
  ntips<-Ntip(tree)
  nedges<-Nedge(tree)
  nodes<-get('nodes',envir=environment(recon.fun),inherits=FALSE)
  ord<-order(nodes[,2])

  disc.recon<-t(Reduce('+',asplit(PP,1))*dx)
  colnames(disc.recon)<-key
  disc.recon<-disc.recon[ord,,drop=FALSE]
  rownames(disc.recon)<-c(tree[['tip.label']],seq_len(tree[['Nnode']])+ntips)

  means<-vars<-matrix(nrow=nedges+1,ncol=k+1)
  dims<-dim(PP)
  for(i in seq_len(k)){
    tmp.P<-matrix(PP[,i,,drop=FALSE],dims[1],dims[3])
    sums<-dx*.colSums(tmp.P,dims[1],dims[3])
    tmp.P<-sweep(tmp.P,2,sums,'/')
    tmp.P[,!sums]<-NA
    means[,i]<-dx*.colSums(sweep(tmp.P,1,xpts,'*'),dims[1],dims[3])
    vars[,i]<-dx*.colSums(sweep(tmp.P,1,xpts.sq,'*'),dims[1],dims[3])-(means[,i])^2
  }
  tmp.P<-Reduce('+',asplit(PP,2))
  sums<-dx*.colSums(tmp.P,dims[1],dims[3])
  tmp.P<-sweep(tmp.P,2,sums,'/')
  means[,k+1]<-dx*.colSums(sweep(tmp.P,1,xpts,'*'),dims[1],dims[3])
  vars[,k+1]<-dx*.colSums(sweep(tmp.P,1,xpts.sq,'*'),dims[1],dims[3])-means[,k+1]^2
  means<-means[ord,,drop=FALSE]
  vars<-vars[ord,,drop=FALSE]
  rownames(means)<-rownames(vars)<-rownames(disc.recon)
  colnames(means)<-colnames(vars)<-c(key,'marginal')
  PP<-PP[,,ord,drop=FALSE]
  dimnames(PP)<-list(NULL,
                     key,
                     rownames(means))
  attr(PP,'xpts')<-NULL
  list('disc.recon'=disc.recon,
       'anc.recon'=list('means'=means,
                        'vars'=vars),
       'ancestral.pdfs'=list('y'=PP,
                             'x'=xpts))
}


