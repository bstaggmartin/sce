#hidden just turns off discrete trait info--make more complicated
#some numerical instability indicates that it might be a good idea to make a "careful" version that directly uses matrix exponentiation
#'@export
make.sce<-function(tree,disc,cont,
                   Q.template=c('ARD','SYM','ER'),
                   exfac=0.5,nx=1024,
                   tip.sig=0,
                   integrate.root=TRUE){
  ##tree topology info##
  list2env(.get.tree.topo.info(tree,pruning=TRUE),envir=environment())
  #continuous
  cont<-cont[tips.ord]
  #discrete
  disc<-disc[tips.ord]
  nas<-is.na(disc)
  disc<-strsplit(disc,'&')
  key<-sort(unique(unlist(disc,use.names=FALSE)))
  k<-length(key)
  disc<-lapply(disc,function(ii) match(ii,key))
  disc[nas]<-list(seq_len(k))

  #Q template
  if(!is.matrix(Q.template)){
    code<-pmatch(Q.template[1],c('ARD','SYM','ER'),nomatch=1)
    Q.template<-matrix(nrow=k,ncol=k)
    if(code==1){
      Q.template[]<-seq_len(k^2)
    }else if(code==2){
      Q.template[lower.tri(Q.template,diag=TRUE)]<-seq_len(k^2/2+k)
      Q.template[upper.tri(Q.template)]<-t(Q.template)[upper.tri(Q.template)]
    }else{
      Q.template[]<-2
      diag(Q.template)<-c(1,3:(k+1))
    }
  }
  nms<-rownames(Q.template)
  if(!is.null(nms)){
    colnames(Q.template)<-nms
    Q.template<-Q.template[key,key]
  }
  par.vec<-sort(unique(as.vector(Q.template)))
  par.vec<-par.vec[par.vec!=0&!is.na(par.vec)]
  Q.template[]<-match(Q.template,par.vec,nomatch=0L)
  par.vec<-rep(0,k^2)
  nz.inds<-Q.template!=0


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
  DFT<-function(x) FFT(c(x,pad),plan=plan)
  #running into some annoying numerical issues--this function might be an inexpensive means of "clean up"
  IDFT<-function(x){
    # out<-Re(IFFT(x,plan=plan))
    # #smarter solution?
    # #set to average of neighbors?
    # #geometric or arithmetic average?
    # #is 0 an appropriate tolerance cutoff?
    # zeros<-out[nonpad]<=0
    # if(any(zeros)){
    #   zeros<-which(zeros)
    #   out[zeros]<-0
    #   l<-out[(zeros-2)%%(2*nx)+1]
    #   r<-out[zeros%%(2*nx)+1]
    #   inds<-l>0&r>0
    #   out[zeros[inds]]<-exp((log(l[inds])+log(r[inds]))/2)
    #   # out[zeros]<-(pmax(0,out[(zeros-2)%%(2*nx)+1])+pmax(0,out[zeros%%(2*nx)+1]))/2
    # }
    # out[nonpad]
    #curious about this...
    #I think this might actually be the way to go!
    abs(IFFT(x,plan=plan)[nonpad])
  }
  #half-"basis" of BM DFT (ignoring other half since it's symmetric)
  basis<- do.call(rbind,rep(list(-0.5*seq(0,pi/dx,length.out=nx+1)^2),k))
  #indices to reflect the above properly
  rep.seq<-c(seq_len(nx+1),nx:2)

  ##formatting trait data into joint array##
  #make array --> dim1 = edges, dim2 = discrete trait, dim3 = continuous trait
  #each cell represents partial conditional likelihood of observed data given that the...
  #...descendant node of a particular edge exhibits a particular discrete-continuous phenotype
  #NEW: try keeping it on frequency spectrum this time?
  #XX is an identically-formatted array used for holding intermediate calculations
  XX<-X<-array(0,c(length(des),k,2*nx))
  scalar.init<- -length(tips)*log(dx)
  #not dealing with hidden stuff for now
  #may want to add support for intraspecific stuff in the future, though...
  new.cont<-findInterval(cont,xpts)
  new.cont<-xpts[new.cont+round((cont-xpts[new.cont])/dx)]
  if(tip.sig){
    tip.sig2<-tip.sig^2
    const<-dnorm(0,0,tip.sig)
    scalar.init<-scalar.init+length(tips)*log(const)
    int.foo<-gauss.DFT(nx,dx,x0=xpts[1])
    foo<-function(x){
      int.foo(x=x,sig2=tip.sig2)/const
    }
  }else{
    foo<-dirac.DFT(nx,dx,x0=xpts[1])
  }
  for(i in seq_along(tips)){
    X[tips[i],disc[[i]],]<-rep(foo(new.cont[i]),each=length(disc[[i]]))
  }

  ##miscellaneous things##
  holder3D<-array(0,c(k,k,nx+1))
  holder2D<-matrix(0,k,nx)
  diag.inds<-rep(c(TRUE,rep(FALSE,k)),length.out=k^2)
  odiag.inds<-!diag.inds
  #about as fast as can be --> just adds ~20ish microseconds for pasting and type conversion
  fast.colsum<-function(x){
    eval(str2lang(paste0('x[',seq_len(dim(x)[1]),',,,drop=FALSE]',collapse='+')))
  }
  fast.rowsum<-function(x){
    eval(str2lang(paste0('x[,',seq_len(dim(x)[2]),',,drop=FALSE]',collapse='+')))
  }

  ##helper functions##
  #converts parameters to convolution function of t
  get.conv.fxn<-function(par){
    par.vec[nz.inds]<-exp(par)[Q.template]
    holder3D[odiag.inds]<-par.vec[odiag.inds]
    holder3D[diag.inds]<- -.rowSums(holder3D[,,1,drop=FALSE],k,k)+par.vec[diag.inds]*basis
    eigs<-apply(holder3D,3,eigen)
    eigvecs<-lapply(eigs,'[[','vectors')
    inveigvecs<-lapply(eigvecs,solve)
    holder3D[]<-unlist(eigvecs,use.names=FALSE)
    eigvecs<-holder3D
    holder3D[]<-unlist(inveigvecs,use.names=FALSE)
    inveigvecs<-holder3D
    holder3D[,1,]<-unlist(lapply(eigs,'[[','values'),use.names=FALSE)
    eigvals<-holder3D[,1,]
    function(t){
      tmp<-sweep(eigvecs,c(2,3),exp(t*eigvals),'*',check.margin=FALSE)
      out<-tmp[,rep.int(1,k),,drop=FALSE]*inveigvecs[rep.int(1,k),,,drop=FALSE]
      for(i in seq_len(k-1)){
        out<-out+tmp[,rep.int(i+1,k),,drop=FALSE]*inveigvecs[rep.int(i+1,k),,,drop=FALSE]
      }
      out[,,rep.seq,drop=FALSE]
    }
  }
  prop.liks<-function(x,t,conv,forwards=FALSE){
    tmp<-conv(t)
    if(forwards){
      fast.colsum(sweep(tmp,c(1,3),x,'*',check.margin=FALSE))
    }else{
      fast.rowsum(sweep(tmp,c(2,3),x,'*',check.margin=FALSE))
    }
  }
  comb.liks<-function(x,G=NULL){
    x<-apply(x,c(1,2),IDFT,simplify=FALSE)
    for(i in seq_len(ncol(x))){
      holder2D[i,]<-Reduce('*',x[,i])
    }
    if(!is.null(G)){
      holder2D<-holder2D*G
    }
    # zeros<-holder2D<1e-10
    # if(any(zeros)){
    #   inds<-which(zeros)
    #   golder2D[zeros]<-(pmin(0,(zeros-2)%%(nx)+1)+pmin(0,zeros%%(nx)+1))/2
    # }
    L<-max(holder2D)
    #check to make sure node has defined conditional likelihoods
    #otherwise some sig2 is too low for accurate calculations...
    if(is.nan(L)|is.infinite(L)|is.na(L)){
      NULL
    }else if(L==0){
      NULL
    }else{
      #rescale conditional likelihoods to prevent underflow
      #also re-DFT
      list(t(apply(holder2D/L,1,DFT)),log(L))
    }
  }
  if(integrate.root){
    root.calc<-function(x,scalar){
      #multiply by scalar, use log-sum-exp trick to sum the conditional likelihoods, weighted by probs they gave rise to data
      #equivalent to Fitzjohn root prior
      #log-sum-exp implicit due to repeated rescaling at nodes...
      #calculation method has been verified to match log(sum(lik^2/sum(lik)))+log(scalar)
      #also, the dx*sum(lik) in the denominator cancels due to multiplication of the entire thing by dx
      x<-apply(x,1,IDFT)
      log(sum(x^2))-log(sum(x))+scalar
    }
  }else{
    root.calc<-function(x,scalar){
      log(max(apply(x,1,IDFT)))+scalar
    }
  }

  #actual function to output
  lik.fun<-function(par,
                    return.array=FALSE){
    conv<-get.conv.fxn(par)
    #scalar --> keeps track of how cond likelihoods are scaled to prevent underflow
    scalar<-scalar.init
    for(i in prune.seq){
      for(j in des[[i]]){
        XX[j,,]<-prop.liks(matrix(X[j,,,drop=FALSE],k,2*nx),elen[j],conv)
      }
      tmp<-comb.liks(XX[des[[i]],,,drop=FALSE])
      if(is.null(tmp)){
        return(Inf)
      }
      X[i,,]<-tmp[[1]]
      scalar<-scalar+tmp[[2]]
    }
    -root.calc(tmp[[1]],scalar)
  }

  recon.fun<-function(par){
    conv<-get.conv.fxn(par)
    for(i in prune.seq){
      for(j in des[[i]]){
        XX[j,,]<-prop.liks(matrix(X[j,,,drop=FALSE],k,2*nx),elen[j],conv)
      }
      tmp<-comb.liks(XX[des[[i]],,,drop=FALSE])
      if(is.null(tmp)){
        stop("Parameters result in undefinied likelihood")
      }
      X[i,,]<-tmp[[1]]
    }
    #X = conditional likelihoods for each edge at TERMINAL NODE given descendants (FFT'd)
    #XX = conditional likelihoods for each edge at STARTING NODE given descendants (FFT'd)
    #G = conditional likelihoods for each edge at TERMINAL NODE given non-descendants (FFT'd)
    G<-X
    G[i,,]<-0;G[i,,1]<-2*nx #root initialization--corresponds to DFT of 1's repped 2*nx times
    #(any phenotype is possible given no information, though may want to allow root priors later on)
    for(i in rev(prune.seq)){
      tmp.G<-t(apply(G[i,,],1,IDFT))
      for(j in des[[i]]){
        tmp<-comb.liks(XX[sis[[j]],,,drop=FALSE],tmp.G)
        if(is.null(tmp)){
          stop("Parameters result in undefined likelihood")
        }
        G[j,,]<-prop.liks(tmp[[1]],elen[j],conv,forwards=TRUE)
      }
    }
    out<-apply(X,c(1,2),IDFT)*apply(G,c(1,2),IDFT)
    #divided by dx to integrate to 1--I think this is right?
    sums<-Reduce('+',asplit(out,c(1,3)))*dx
    out<-aperm(sweep(out,2,sums,'/'),c(1,3,2))
    attr(out,'xpts')<-xpts
    out
  }

  out<-list('lik.fun'=lik.fun,'recon.fun'=recon.fun)
  attr(out,'key')<-key
  attr(out,'Q.template')<-Q.template
  out
}

#' #hidden just turns off discrete trait info--make more complicated
#' #'@export
#' make.sce<-function(tree,disc,cont,
#'                    exfac=0.5,nx=1024,nt=100,
#'                    tip.sig=0,
#'                    integrate.root=TRUE,
#'                    hidden=FALSE){
#'   ##tree topology info##
#'   #reorder for pruning
#'   attr(tree,'order')<-NULL
#'   tree<-reorder(tree,order='pruningwise')
#'   #get ancestors
#'   anc<-match(tree$edge[,1],tree$edge[,2])
#'   nedge<-length(anc)+1
#'   anc[is.na(anc)]<-nedge #root edge index is number of edges + 1
#'   #get descendants, note where tips are, make pruning sequence
#'   sis<-des<-rep(list(integer(0)),nedge)
#'   prune.seq<-seq_len(nedge)
#'   tmp<-split(prune.seq[-nedge],anc)
#'   des[as.numeric(names(tmp))]<-tmp
#'   tips<-which(lengths(des)==0)
#'   prune.seq<-prune.seq[-tips]
#'   #sister edges...more complicated, but needed for ancestral state recon
#'   ndes<-lengths(des)
#'   di<-ndes==2
#'   di.des<-unlist(des[di])
#'   if(length(di.des)){
#'     odds<-seq.int(1,length(di.des),2)
#'     evens<-odds+1
#'     odds<-di.des[odds]
#'     evens<-di.des[evens]
#'     sis[odds]<-evens
#'     sis[evens]<-odds
#'   }
#'   poly.des<-des[!di]
#'   if(length(poly.des)){
#'     ndes<-ndes[!di]
#'     unlist.poly.des<-unlist(poly.des,use.names=FALSE)
#'     sis[unlist.poly.des]<-rep(poly.des,ndes)
#'     foo<-function(x){
#'       tmp<-sis[[x]]
#'       tmp[tmp!=x]
#'     }
#'     sis[unlist.poly.des]<-lapply(unlist.poly.des,foo)
#'   }
#'   #reorder tip data
#'   tips.ord<-tree$tip.label[tree$edge[tips,2]]
#'   #continuous
#'   cont<-cont[tips.ord]
#'   #discrete
#'   disc<-disc[tips.ord]
#'   if(!is.factor(disc)){
#'     disc<-as.factor(disc)
#'   }
#'   key<-levels(disc)
#'   k<-length(key)
#'   disc<-as.numeric(disc)
#'   if(hidden){
#'     disc<-rep(1,length(tips))
#'     key<-letters[seq_len(k)]
#'   }
#'
#'   ##discretizing time##
#'   #get tree height and edge lengths
#'   hgt<-max(node.depth.edgelength(tree))
#'   elen<-c(tree$edge.length,0) #extra 0 for root edge
#'   #"ideal" time interval based on nt
#'   dt<-hgt/nt
#'   #find time res along each branch
#'   nts<-ceiling(elen/dt)
#'   #actual time interval along each branch
#'   dts<-elen/nts
#'   dts[is.nan(dts)]<-0 #for 0-length branches (like root edge)
#'   nts<-lapply(nts,seq_len)
#'
#'   ##discretizing continuous trait##
#'   #make a DFT plan
#'   plan<-planFFT(2*nx)
#'   #find min and max considered values for continuous trait
#'   xlim<-range(cont)+c(-1,1)*exfac*diff(range(cont))
#'   #get the breaks and find the inter-break distance
#'   xpts<-seq(xlim[1],xlim[2],length.out=nx)
#'   dx<-diff(xpts[c(1,2)])
#'   #make a vector of 0s for padding
#'   pad<-rep(0,nx)
#'   nonpad<-rep(c(TRUE,FALSE),each=nx)
#'   #function for getting/inverting DFTs of conditional likelihoods
#'   DFT.X<-function(x) FFT(c(x,pad),plan=plan)
#'   IDFT.X<-function(X) Re(apply(X,2,IFFT,plan=plan)[nonpad,,drop=FALSE])
#'   #make matrix of DFTs for normal distributions
#'   gauss.DFT<-function(nx,dx){
#'     tmp<-c(0,seq_len(nx))
#'     out<-c((-2*pi^2*(tmp/(2*nx))^2),(-2*pi^2*((tmp[-c(1,2)]+nx-1)/(2*nx)-1)^2))
#'     out/dx^2
#'   }
#'   norms<-matrix(gauss.DFT(nx,dx),nrow=2*nx,ncol=k)
#'   #function for scaling normal distribution DFTs based on sig2/time interval
#'   DFT.norms<-function(scalar) exp(sweep(norms,2,scalar,'*',check.margin=FALSE))
#'
#'   ##formatting trait data into joint array##
#'   #make array --> dim1 = edges, dim2 = continuous trait, dim3 = discrete trait
#'   #each cell represents partial conditional likelihood of observed data given that the...
#'   #...descendant node of a particular edge exhibits a particular discrete-continuous phenotype
#'   X<-array(0,c(length(des),nx,k))
#'   #XX is an identically-formatted array used for holding intermediate calculations
#'   XX<-X
#'   scalar.init<- -length(tips)*log(dx)
#'   if(!tip.sig){
#'     inds<-findInterval(cont,xpts)
#'     inds<-inds+round((cont-xpts[inds])/dx)
#'     inds<-cbind(tips,inds,disc)
#'     X[inds]<-1
#'     if(hidden){
#'       X[]<-X[,,1,drop=FALSE]
#'     }
#'   }else{
#'     for(i in seq_along(tips)){
#'       tmp<-dnorm(xpts,cont[i],tip.sig)
#'       tmp<-tmp/sum(tmp) #helps prevent aliasing issues
#'       tmp.max<-max(tmp)
#'       scalar.init<-scalar.init+log(tmp.max)
#'       tmp<-tmp/tmp.max
#'       if(hidden){
#'         X[tips[i],,]<-tmp
#'       }else{
#'         X[tips[i],,disc[i]]<-tmp
#'       }
#'     }
#'   }
#'   if(hidden){
#'     scalar.init<-scalar.init+(k-1)*scalar.init #add extra scalars for extra conditional likelihoods...
#'   }
#'   scalar.init<-scalar.init-1 #no idea where extra -1 comes from
#'
#'   ##miscellaneous things##
#'   QQ<-Q<-matrix(0,nrow=k,ncol=k)
#'   Q.diag<-seq.int(1,k^2,k+1)
#'   Q.inds<-seq_len(k^2)[-Q.diag]
#'   sig2.par<-rep(TRUE,k^2)
#'   sig2.par[seq_len(k^2-k)]<-FALSE
#'   #about as fast as can be --> just adds ~20ish microseconds for pasting and type conversion
#'   mult<-function(XX){
#'     eval(str2lang(paste0('XX[',seq_len(dim(XX)[1]),',,,drop=FALSE]',collapse='*')))
#'   }
#'
#'   ##helper functions##
#'   form.par<-function(par,forwards=FALSE){
#'     #format parameters --> given on log scale since they all must be positive
#'     par<-exp(par)
#'     Q[Q.inds]<-par[!sig2.par]
#'     Q[Q.diag]<- -.rowSums(Q,k,k)
#'     sig2<-par[sig2.par]
#'     #eigen decomp of Q
#'     eig<-eigen(Q)
#'     eigval<-eig$values
#'     eigvec<-eig$vectors
#'     inveigvec<-solve(eigvec)
#'     int.get.P<-function(t){
#'       #transition probability matrix for j's time interval --> expm(Q*dts[j])
#'       #if Ve = eigenvectors of Q and Va = eigenvalues, then Ve*diag(exp(dts[j]*Va))*Ve^-1
#'       #transpose prob mat because we're calculating conditional likelihoods, not probs
#'       QQ[Q.diag]<-exp(eigval*t)
#'       Re(eigvec%*%QQ%*%inveigvec)
#'     }
#'     #if this is backwards, P should be backwards (transpose it to do this)
#'     if(forwards){
#'       get.P<-int.get.P
#'     }else{
#'       get.P<-function(t){
#'         t(int.get.P(t))
#'       }
#'     }
#'     list(sig2,get.P)
#'   }
#'   prop.liks<-function(X,get.P,sig2,t,nt,forwards=FALSE){
#'     #get DFTs of conditional likelihoods
#'     X<-apply(X,2,DFT.X)
#'     #get probability transition matrix for time interval
#'     P<-get.P(t)
#'     #scale normal DFTs for each state given sig2 and time interval
#'     N<-DFT.norms(sig2*t)
#'     if(forwards){
#'       for(t in nt){
#'         X<-X*N
#'         X<-X%*%P
#'       }
#'     }else{
#'       for(t in nt){
#'         #use transition probability matrix to shuffle discrete states around
#'         X<-X%*%P
#'         #use normal DFTs for convolution, shuffling continuous states around
#'         X<-X*N
#'       }
#'     }
#'     #invert DFT
#'     IDFT.X(X)
#'   }
#'   comb.liks<-function(X){
#'     #mutliply conditional likelihoods to get conditional likelihoods at focal node
#'     X<-mult(X)
#'     L<-max(X)
#'     #check to make sure node has defined conditional likelihoods
#'     #otherwise some sig2 is too low for accurate calculations...
#'     if(is.nan(L)|is.infinite(L)|is.na(L)){
#'       NULL
#'     }else if(L<0){
#'       NULL
#'     }else{
#'       #rescale conditional likelihoods to prevent underflow
#'       list(X/L,log(L))
#'     }
#'   }
#'   if(integrate.root){
#'     root.calc<-function(lik,scalar){
#'       #multiply by scalar, use log-sum-exp trick to sum the conditional likelihoods, weighted by probs they gave rise to data
#'       #equivalent to Fitzjohn root prior
#'       #log-sum-exp implicit due to repeated rescaling at nodes...
#'       #calculation method has been verified to match log(sum(lik^2/sum(lik)))+log(scalar)
#'       #also, the dx*sum(lik) in the denominator cancels due to multiplication of the entire thing by dx
#'       lik[lik<0]<-0
#'       log(sum(lik^2))-log(sum(lik))+scalar
#'     }
#'   }else{
#'     root.calc<-function(lik,scalar){
#'       max(lik)+scalar
#'     }
#'   }
#'
#'   #actual function to output
#'   lik.fun<-function(par,
#'                     return.array=FALSE){
#'     tmp<-form.par(par)
#'     sig2<-tmp[[1]]
#'     get.P<-tmp[[2]]
#'     #scalar --> keeps track of how cond likelihoods are scaled to prevent underflow
#'     scalar<-scalar.init
#'     for(i in prune.seq){
#'       tmp.des<-des[[i]]
#'       for(j in tmp.des){
#'         XX[j,,]<-prop.liks(matrix(X[j,,,drop=FALSE],nx,k),get.P,sig2,dts[j],nts[[j]])
#'       }
#'       tmp<-comb.liks(XX[tmp.des,,,drop=FALSE])
#'       if(is.null(tmp)){
#'         return(Inf)
#'       }
#'       X[i,,]<-tmp[[1]]
#'       scalar<-scalar+tmp[[2]]
#'     }
#'     if(return.array){
#'       list(X,XX)
#'     }else{
#'       -root.calc(tmp[[1]],scalar)
#'     }
#'   }
#'
#'   recon.fun<-function(par){
#'     tmp<-form.par(par,forwards=TRUE)
#'     sig2<-tmp[[1]]
#'     get.P<-tmp[[2]]
#'     tmp<-lik.fun(par,return.array=TRUE)
#'     #X and XX = conditional likelihoods for each edge given descendants
#'     #XX contains conditional likelihood propagated down edge towards root
#'     #G = conditional likelihoods for each edge given non-descendants
#'     X<-tmp[[1]]
#'     XX<-tmp[[2]]
#'     G<-X
#'     G[nedge,,]<-1 #root initialization
#'     for(i in seq.int(nedge-1,1)){
#'       tmp.anc<-anc[i]
#'       tmp.sis<-sis[[i]]
#'       tmp.G<-mult(XX[tmp.sis,,,drop=FALSE])
#'       tmp.G<-tmp.G*G[tmp.anc,,,drop=FALSE]
#'       tmp.max<-max(tmp.G)
#'       tmp.G<-tmp.G/tmp.max
#'       G[i,,]<-prop.liks(matrix(tmp.G,nx,k),get.P,sig2,dts[i],nts[[i]],
#'                         forwards=TRUE)
#'     }
#'     PP<-G*X
#'     #divided by dx to integrate to 1
#'     sums<-Reduce('+',asplit(PP,c(2,3)))*dx
#'     PP<-sweep(PP,1,sums,'/')
#'     attr(PP,'xpts')<-xpts
#'     PP
#'   }
#'
#'   out<-list('lik.fun'=lik.fun,'recon.fun'=recon.fun)
#'   attr(out,'key')<-key
#'   out
#' }
#'
#' #'@export
#' make.sce2<-function(tree,disc,cont,
#'                    exfac=0.5,nx=1024,
#'                    tip.sig=0,
#'                    integrate.root=TRUE){
#'   ##tree topology info##
#'   list2env(.get.tree.topo.info(tree,pruning=TRUE),envir=environment())
#'   #continuous
#'   cont<-cont[tips.ord]
#'   #discrete
#'   disc<-disc[tips.ord]
#'   if(!is.factor(disc)){
#'     disc<-as.factor(disc)
#'   }
#'   key<-levels(disc)
#'   k<-length(key)
#'   disc<-as.numeric(disc)
#'
#'   ##discretizing continuous trait##
#'   #find min and max considered values for continuous trait
#'   xlim<-range(cont)+c(-1,1)*exfac*(max(cont)-min(cont))
#'   #get the breaks and find the inter-break distance
#'   xpts<-seq(xlim[1],xlim[2],length.out=nx)
#'   dx<-xpts[2]-xpts[1]
#'
#'   ##some misc initialization##
#'   nx.seq<-seq_len(nx)
#'   k.seq<-seq_len(k)
#'   knx<-k*nx
#'   knx.seq<-seq_len(knx)
#'   invdxsq<-1/(2*dx^2)
#'
#'   ##formatting trait data into joint array##
#'   #make array --> dim1 = trait value, dim2 = edges
#'   #each cell represents partial conditional likelihood of observed data given that the...
#'   #...descendant node of a particular edge exhibits a particular discrete-continuous phenotype
#'   X<-matrix(1,knx,length(des))
#'   X[,tips]<-0
#'   scalar.init<- -length(tips)*log(dx)
#'   if(!tip.sig){
#'     inds<-findInterval(cont,xpts)
#'     inds<-inds+round((cont-xpts[inds])/dx)
#'     inds<-cbind(inds+(disc-1)*nx,tips)
#'     X[inds]<-1
#'   }else{
#'     for(i in seq_along(tips)){
#'       tmp<-dnorm(xpts,cont[i],tip.sig)
#'       tmp<-tmp/sum(tmp) #helps prevent aliasing issues
#'       tmp.max<-max(tmp)
#'       scalar.init<-scalar.init+log(tmp.max)
#'       tmp<-tmp/tmp.max
#'       X[nx.seq+(disc[i]-1)*nx,tips[i]]<-tmp
#'     }
#'   }
#'   scalar.init<-scalar.init-1 #no idea where extra -1 comes from
#'
#'   #make convolution matrix fxn
#'   rate.inds<-cbind(nx.seq[-1],nx.seq[-nx])
#'   rate.inds<-rbind(rate.inds,rate.inds[,c(2,1)])
#'   tran.inds<-lapply(k.seq,function(ii) nx.seq+(ii-1)*nx)
#'   inds<-matrix(list(),k,k)
#'   for(j in k.seq){
#'     for(i in k.seq){
#'       if(i==j){
#'         inds[[i,j]]<-rate.inds+(i-1)*nx
#'       }else{
#'         inds[[i,j]]<-do.call(cbind,tran.inds[c(i,j)])
#'       }
#'     }
#'   }
#'   holder<-matrix(0,knx,knx)
#'   diag.inds<-matrix(knx.seq,knx,2)
#'   sig2.inds<-seq.int(1,k^2,k+1)
#'
#'   ##helper functions##
#'   get.Q<-function(par){
#'     par<-exp(par)
#'     par[sig2.inds]<-par[sig2.inds]*invdxsq
#'     for(i in seq_along(par)){
#'       holder[inds[[i]]]<-par[i]
#'     }
#'     holder[diag.inds]<- -.rowSums(holder,knx,knx)
#'     holder
#'   }
#'   if(integrate.root){
#'     root.calc<-function(lik,scalar){
#'       #multiply by scalar, use log-sum-exp trick to sum the conditional likelihoods, weighted by probs they gave rise to data
#'       #equivalent to Fitzjohn root prior
#'       #log-sum-exp implicit due to repeated rescaling at nodes...
#'       #calculation method has been verified to match log(sum(lik^2/sum(lik)))+log(scalar)
#'       #also, the dx*sum(lik) in the denominator cancels due to multiplication of the entire thing by dx
#'       lik[lik<0]<-0
#'       log(sum(lik^2))-log(sum(lik))+scalar
#'     }
#'   }else{
#'     root.calc<-function(lik,scalar){
#'       max(lik)+scalar
#'     }
#'   }
#'
#'   #actual function to output
#'   lik.fun<-function(par){
#'     Q<-get.Q(par)
#'     #scalar --> keeps track of how cond likelihoods are scaled to prevent underflow
#'     scalar<-scalar.init
#'     for(i in prune.seq){
#'       tmp.des<-des[[i]]
#'       for(j in tmp.des){
#'         X[,i]<-X[,i]*expAtv(Q,X[,j],elen[j])[[1]]
#'       }
#'       L<-max(X[,i])
#'       if(is.nan(L)|is.infinite(L)|is.na(L)){
#'         return(Inf)
#'       }else if(L<0){
#'         return(Inf)
#'       }
#'       X[,i]<-X[,i]/L
#'       scalar<-scalar+log(L)
#'     }
#'     -root.calc(X[,i],scalar)
#'   }
#'
#'   out<-list('lik.fun'=lik.fun)
#'   attr(out,'key')<-key
#'   out
#' }
