#hidden just turns off discrete trait info--make more complicated
#'@export
make.sce<-function(tree,disc,cont,
                   exfac=0.5,nx=1024,nt=100,
                   tip.sig=0,
                   integrate.root=TRUE,
                   hidden=FALSE){
  ##tree topology info##
  #reorder for pruning
  attr(tree,'order')<-NULL
  tree<-reorder(tree,order='pruningwise')
  #get ancestors
  anc<-match(tree$edge[,1],tree$edge[,2])
  nedge<-length(anc)+1
  anc[is.na(anc)]<-nedge #root edge index is number of edges + 1
  #get descendants, note where tips are, make pruning sequence
  sis<-des<-rep(list(integer(0)),nedge)
  prune.seq<-seq_len(nedge)
  tmp<-split(prune.seq[-nedge],anc)
  des[as.numeric(names(tmp))]<-tmp
  tips<-which(lengths(des)==0)
  prune.seq<-prune.seq[-tips]
  #sister edges...more complicated, but needed for ancestral state recon
  ndes<-lengths(des)
  di<-ndes==2
  di.des<-unlist(des[di])
  if(length(di.des)){
    odds<-seq.int(1,length(di.des),2)
    evens<-odds+1
    odds<-di.des[odds]
    evens<-di.des[evens]
    sis[odds]<-evens
    sis[evens]<-odds
  }
  poly.des<-des[!di]
  if(length(poly.des)){
    ndes<-ndes[!di]
    unlist.poly.des<-unlist(poly.des,use.names=FALSE)
    sis[unlist.poly.des]<-rep(poly.des,ndes)
    foo<-function(x){
      tmp<-sis[[x]]
      tmp[tmp!=x]
    }
    sis[unlist.poly.des]<-lapply(unlist.poly.des,foo)
  }
  #reorder tip data
  tips.ord<-tree$tip.label[tree$edge[tips,2]]
  #continuous
  cont<-cont[tips.ord]
  #discrete
  disc<-disc[tips.ord]
  if(!is.factor(disc)){
    disc<-as.factor(disc)
  }
  key<-levels(disc)
  k<-length(key)
  disc<-as.numeric(disc)
  if(hidden){
    disc<-rep(1,n)
    key<-letters[seq_len(k)]
  }
  
  ##discretizing time##
  #get tree height and edge lengths
  hgt<-max(node.depth.edgelength(tree))
  elen<-c(tree$edge.length,0) #extra 0 for root edge
  #"ideal" time interval based on nt
  dt<-hgt/nt
  #find time res along each branch
  nts<-ceiling(elen/dt)
  #actual time interval along each branch
  dts<-elen/nts
  dts[is.nan(dts)]<-0 #for 0-length branches (like root edge)
  nts<-lapply(nts,seq_len)
  
  ##discretizing continuous trait##
  #make a DFT plan
  plan<-planFFT(2*nx)
  #find min and max considered values for continuous trait
  xlim<-range(cont)+c(-1,1)*exfac*diff(range(cont))
  #get the breaks and find the inter-break distance
  xpts<-seq(xlim[1],xlim[2],length.out=nx)
  dx<-diff(xpts[c(1,2)])
  #make a vector of 0s for padding
  pad<-rep(0,nx)
  nonpad<-rep(c(TRUE,FALSE),each=nx)
  #function for getting/inverting DFTs of conditional likelihoods
  DFT.X<-function(x) FFT(c(x,pad),plan=plan)
  IDFT.X<-function(X) Re(apply(X,2,IFFT,plan=plan)[nonpad,,drop=FALSE])
  #make matrix of DFTs for normal distributions
  gauss.DFT<-function(nx,dx){
    tmp<-c(0,seq_len(nx))
    out<-c((-2*pi^2*(tmp/(2*nx))^2),(-2*pi^2*((tmp[-c(1,2)]+nx-1)/(2*nx)-1)^2))
    out/dx^2
  }
  norms<-matrix(gauss.DFT(nx,dx),nrow=2*nx,ncol=k)
  #function for scaling normal distribution DFTs based on sig2/time interval
  DFT.norms<-function(scalar) exp(sweep(norms,2,scalar,'*',check.margin=FALSE))
  
  ##formatting trait data into joint array##
  #make array --> dim1 = edges, dim2 = continuous trait, dim3 = discrete trait
  #each cell represents partial conditional likelihood of observed data given that the...
  #...descendant node of a particular edge exhibits a particular discrete-continuous phenotype
  X<-array(0,c(length(des),nx,k))
  #XX is an identically-formatted array used for holding intermediate calculations
  XX<-X
  scalar.init<- -length(tips)*log(dx)
  if(!tip.sig){
    inds<-findInterval(cont,xpts)
    inds<-inds+round((cont-xpts[inds])/dx)
    inds<-cbind(tips,inds,disc)
    X[inds]<-1
    if(hidden){
      X[]<-X[,,1,drop=FALSE]
    }
  }else{
    for(i in seq_along(tips)){
      tmp<-dnorm(xpts,cont[i],tip.sig)
      tmp<-tmp/sum(tmp) #helps prevent aliasing issues
      tmp.max<-max(tmp)
      scalar.init<-scalar.init+log(tmp.max)
      tmp<-tmp/tmp.max
      if(hidden){
        X[tips[i],,]<-tmp
      }else{
        X[tips[i],,disc[i]]<-tmp
      }
    }
  }
  if(hidden){
    scalar.init<-scalar.init+(k-1)*scalar.init #add extra scalars for extra conditional likelihoods...
  }
  scalar.init<-scalar.init-1 #no idea where extra -1 comes from
  
  ##miscellaneous things##
  QQ<-Q<-matrix(0,nrow=k,ncol=k)
  Q.diag<-seq.int(1,k^2,k+1)
  Q.inds<-seq_len(k^2)[-Q.diag]
  sig2.par<-rep(TRUE,k^2)
  sig2.par[seq_len(k^2-k)]<-FALSE
  #about as fast as can be --> just adds ~20ish microseconds for pasting and type conversion
  mult<-function(XX){
    eval(str2lang(paste0('XX[',seq_len(dim(XX)[1]),',,,drop=FALSE]',collapse='*')))
  }
  
  ##helper functions##
  form.par<-function(par,forwards=FALSE){
    #format parameters --> given on log scale since they all must be positive
    par<-exp(par)
    Q[Q.inds]<-par[!sig2.par]
    Q[Q.diag]<- -.rowSums(Q,k,k)
    sig2<-par[sig2.par]
    #eigen decomp of Q
    eig<-eigen(Q)
    eigval<-eig$values
    eigvec<-eig$vectors
    inveigvec<-solve(eigvec)
    int.get.P<-function(t){
      #transition probability matrix for j's time interval --> expm(Q*dts[j])
      #if Ve = eigenvectors of Q and Va = eigenvalues, then Ve*diag(exp(dts[j]*Va))*Ve^-1
      #transpose prob mat because we're calculating conditional likelihoods, not probs
      QQ[Q.diag]<-exp(eigval*t)
      Re(eigvec%*%QQ%*%inveigvec)
    }
    #if this is backwards, P should be backwards (transpose it to do this)
    if(forwards){
      get.P<-int.get.P
    }else{
      get.P<-function(t){
        t(int.get.P(t))
      }
    }
    list(sig2,get.P)
  }
  prop.liks<-function(X,get.P,sig2,t,nt,forwards=FALSE){
    #get DFTs of conditional likelihoods
    X<-apply(X,2,DFT.X)
    #get probability transition matrix for time interval
    P<-get.P(t)
    #scale normal DFTs for each state given sig2 and time interval
    N<-DFT.norms(sig2*t)
    if(forwards){
      for(t in nt){
        X<-X*N
        X<-X%*%P
      }
    }else{
      for(t in nt){
        #use transition probability matrix to shuffle discrete states around
        X<-X%*%P
        #use normal DFTs for convolution, shuffling continuous states around
        X<-X*N
      }
    }
    #invert DFT
    IDFT.X(X)
  }
  comb.liks<-function(X){
    #mutliply conditional likelihoods to get conditional likelihoods at focal node
    X<-mult(X)
    L<-max(X)
    #check to make sure node has defined conditional likelihoods
    #otherwise some sig2 is too low for accurate calculations...
    if(is.nan(L)|is.infinite(L)|is.na(L)){
      NULL
    }else if(L<0){
      NULL
    }else{
      #rescale conditional likelihoods to prevent underflow
      list(X/L,log(L))
    }
  }
  if(integrate.root){
    root.calc<-function(lik,scalar){
      #multiply by scalar, use log-sum-exp trick to sum the conditional likelihoods, weighted by probs they gave rise to data
      #equivalent to Fitzjohn root prior
      #log-sum-exp implicit due to repeated rescaling at nodes...
      #calculation method has been verified to match log(sum(lik^2/sum(lik)))+log(scalar)
      #also, the dx*sum(lik) in the denominator cancels due to multiplication of the entire thing by dx
      lik[lik<0]<-0
      log(sum(lik^2))-log(sum(lik))+scalar
    }
  }else{
    root.calc<-function(lik,scalar){
      max(lik)+scalar
    }
  }
  
  #actual function to output
  lik.fun<-function(par,
                    return.array=FALSE){
    tmp<-form.par(par)
    sig2<-tmp[[1]]
    get.P<-tmp[[2]]
    #scalar --> keeps track of how cond likelihoods are scaled to prevent underflow
    scalar<-scalar.init
    for(i in prune.seq){
      tmp.des<-des[[i]]
      for(j in tmp.des){
        XX[j,,]<-prop.liks(matrix(X[j,,,drop=FALSE],nx,k),get.P,sig2,dts[j],nts[[j]])
      }
      tmp<-comb.liks(XX[tmp.des,,,drop=FALSE])
      if(is.null(tmp)){
        return(Inf)
      }
      X[i,,]<-tmp[[1]]
      scalar<-scalar+tmp[[2]]
    }
    if(return.array){
      list(X,XX)
    }else{
      -root.calc(tmp[[1]],scalar)
    }
  }
  
  recon.fun<-function(par){
    tmp<-form.par(par,forwards=TRUE)
    sig2<-tmp[[1]]
    get.P<-tmp[[2]]
    tmp<-lik.fun(par,return.array=TRUE)
    #X and XX = conditional likelihoods for each edge given descendants
    #XX contains conditional likelihood propagated down edge towards root
    #G = conditional likelihoods for each edge given non-descendants
    X<-tmp[[1]]
    XX<-tmp[[2]]
    G<-X
    G[nedge,,]<-1 #root initialization
    for(i in seq.int(nedge-1,1)){
      tmp.anc<-anc[i]
      tmp.sis<-sis[[i]]
      tmp.G<-mult(XX[tmp.sis,,,drop=FALSE])
      tmp.G<-tmp.G*G[tmp.anc,,,drop=FALSE]
      tmp.max<-max(tmp.G)
      tmp.G<-tmp.G/tmp.max
      G[i,,]<-prop.liks(matrix(tmp.G,nx,k),get.P,sig2,dts[i],nts[[i]],
                        forwards=TRUE)
    }
    PP<-G*X
    #divided by dx to integrate to 1
    sums<-Reduce('+',asplit(PP,c(2,3)))*dx
    PP<-sweep(PP,1,sums,'/')
    attr(PP,'xpts')<-xpts
    PP
  }
  
  out<-list('lik.fun'=lik.fun,'recon.fun'=recon.fun)
  attr(out,'key')<-key
  out
}
