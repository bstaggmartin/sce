#temporary maximum likelihood search function
#' @export
find.mle<-function(sce,inits=NULL,...,
                   inits.mean=1,inits.sd=10){
  key<-attr(sce,'key')
  k<-length(key)

  lik.fun<-sce[['lik.fun']]

  args.ls<-list(...)
  if(is.null(args.ls$method)){
    args.ls$method<-'BFGS'
  }
  Q.template<-get('Q.template',envir=environment(lik.fun),inherits=FALSE)
  if(is.null(inits)){
    inits<-rnorm(max(Q.template),log(inits.mean),log(inits.sd))
    # inits<-vector('numeric',k^2)
    # tree<-get('tree',envir=environment(lik.fun),inherits=FALSE)
    # hgt<-max(node.depth.edgelength(tree))
    # xlim<-get('xlim',envir=environment(lik.fun),inherits=FALSE)
    # diag.inds<-get('diag.inds',envir=environment(lik.fun),inherits=FALSE)
    # #initial stab at default inits that prevent underflow
    # #bonus about now using absolute value of inverse fourier transform--underflow is no longer a concern (though numerical precision probably gets wonky...)
    # inits[!diag.inds]<-log(rexp(sum(!diag.inds),1/hgt)) #average 1 transition for lineage evolving from root to tips
    # inits[diag.inds]<-2*rnorm(sum(diag.inds),log(diff(xlim)/4),log(1.5)) #~50% increases/reductions in sd reasonable
  }
  inits<-rep(inits,length.out=max(Q.template))
  args.ls$par<-inits
  args.ls$fn<-lik.fun

  res<-do.call(optim,args.ls)

  diag.inds<-get('diag.inds',envir=environment(lik.fun),inherits=FALSE)
  nz.inds<-get('nz.inds',envir=environment(lik.fun),inherits=FALSE)
  par.vec<-rep(0,k^2)
  par.vec[nz.inds]<-exp(res$par)[Q.template]
  Q<-matrix(0,k,k,dimnames=list(key,key))
  Q[!diag.inds]<-par.vec[!diag.inds]
  diag(Q)<- -rowSums(Q)
  sig2<-setNames(par.vec[diag.inds],key)
  lik<- -res$value

  list('Q'=Q,
       'sig2'=sig2,
       'lnL'=lik,
       'sce'=sce,
       'par'=res$par,
       'hessian'=if(!is.null(res$hessian)) res$hessian else NULL)
}
