#temporary maximum likelihood search function
#' @export
find.mle<-function(sce,inits=NULL,...){
  key<-attr(sce,'key')
  k<-length(key)
  
  lik.fun<-sce[['lik.fun']]
  
  args.ls<-list(...)
  if(is.null(args.ls$method)){
    args.ls$method<-'BFGS'
  }
  if(is.null(inits)){
    inits<-vector('numeric',k^2)
    hgt<-get('hgt',envir=environment(lik.fun),inherits=FALSE)
    xlim<-get('xlim',envir=environment(lik.fun),inherits=FALSE)
    sig2.par<-get('sig2.par',envir=environment(lik.fun),inherits=FALSE)
    #initial stab at default inits that prevent underflow
    inits[!sig2.par]<-log(rexp(k^2-k,1/hgt)) #average 1 transition for lineage evolving from root to tips
    inits[sig2.par]<-2*rnorm(k,log(diff(xlim)/2),log(1.5)) #~1/2 of x range, ~50% increases/reductions in sd reasonable
  }
  inits<-rep(inits,length.out=k^2)
  args.ls$par<-inits
  args.ls$fn<-lik.fun
  
  res<-do.call(optim,args.ls)
  
  par<-exp(res$par)
  Q<-matrix(0,k,k,dimnames=list(key,key))
  Q[seq_len(k^2)[-seq.int(1,k^2,k+1)]]<-par[seq_len(k^2-k)]
  diag(Q)<- -rowSums(Q)
  sig2<-setNames(par[seq_len(k)+k^2-k],key)
  lik<- -res$value
  
  list('Q'=Q,
       'sig2'=sig2,
       'lnL'=lik,
       'sce'=sce,
       'par'=res$par)
}
