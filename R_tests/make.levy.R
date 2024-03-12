#' @export
make.levy<-function(tree,dat,model=c('JN','NIG','VG'),
                    exfac=0.5,nx=1024,nt=100,
                    tip.sig=0,
                    integrate.root=TRUE){
  list2env(.get.tree.topo.info(tree,pruning=TRUE),envir=environment())
  choices<-c('JN','NIG','VG')
  if(is.character(model)){
    model<-pmatch(model[1],choices)
  }
  cont<-dat[tips.ord]

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
  DFT.X<-function(X) FFT(c(X,pad),plan=plan)
  IDFT.X<-function(X) Re(IFFT(X,plan=plan)[nonpad])
  #make DFTs for transition probabilities
  DFT.N<-BM.DFT(nx,dx)
  DFT.L<-switch(model,
                JN.DFT(nx,dx),
                NIG.DFT(nx,dx),
                VG.DFT(nx,dx))

  ##formatting trait data into joint array##
  #make array --> dim1 = edges, dim2 = continuous trait
  #each cell represents partial conditional likelihood of observed data given that the...
  #...descendant node of a particular edge exhibits a particular discrete-continuous phenotype
  X<-matrix(0,length(des),nx)
  #XX is an identically-formatted array used for holding intermediate calculations
  XX<-X
  scalar.init<- -length(tips)*log(dx)
  if(!tip.sig){
    inds<-findInterval(cont,xpts)
    inds<-inds+round((cont-xpts[inds])/dx)
    inds<-cbind(tips,inds)
    X[inds]<-1
  }else{
    for(i in seq_along(tips)){
      tmp<-dnorm(xpts,cont[i],tip.sig)
      tmp<-tmp/sum(tmp) #helps prevent aliasing issues
      tmp.max<-max(tmp)
      scalar.init<-scalar.init+log(tmp.max)
      tmp<-tmp/tmp.max
      X[tips[i],]<-tmp
    }
  }
  scalar.init<-scalar.init-1 #no idea where extra -1 comes from

  ##miscellaneous things##
  #about as fast as can be --> just adds ~20ish microseconds for pasting and type conversion
  mult<-function(XX){
    eval(str2lang(paste0('XX[',seq_len(dim(XX)[1]),',,drop=FALSE]',collapse='*')))
  }

  ##helper functions##
  prop.liks<-function(X,rate,jump,sig2,t){
    #get DFTs of conditional likelihoods
    X<-DFT.X(X)
    #levy DFT for given t
    L<-DFT.L(t,rate,jump)
    #normal DFT for given t
    N<-DFT.N(t,sig2)
    #invert convolved DFT
    IDFT.X(X*L*N)
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
    #par<-exp(par)
    if(any(par<0)){
      return(Inf)
    }
    rate<-par[1]
    jump<-par[2]
    sig2<-par[3]
    #scalar --> keeps track of how cond likelihoods are scaled to prevent underflow
    scalar<-scalar.init
    for(i in prune.seq){
      tmp.des<-des[[i]]
      for(j in tmp.des){
        XX[j,]<-prop.liks(X[j,],rate,jump,sig2,elen[j])
      }
      tmp<-comb.liks(XX[tmp.des,,drop=FALSE])
      if(is.null(tmp)){
        return(Inf)
      }
      X[i,]<-tmp[[1]]
      scalar<-scalar+tmp[[2]]
    }
    if(return.array){
      list(X,XX)
    }else{
      -root.calc(tmp[[1]],scalar)
    }
  }

  recon.fun<-function(par){
    tmp<-lik.fun(par,return.array=TRUE)
    if(!is.list(tmp)){
      stop('underflow')
    }
    #par<-exp(par)
    rate<-par[1]
    jump<-par[2]
    sig2<-par[3]
    #X and XX = conditional likelihoods for each edge given descendants
    #XX contains conditional likelihood propagated down edge towards root
    #G = conditional likelihoods for each edge given non-descendants
    X<-tmp[[1]]
    XX<-tmp[[2]]
    G<-X
    G[nedge,]<-1 #root initialization
    for(i in seq.int(nedge-1,1)){
      tmp.anc<-anc[i]
      tmp.sis<-sis[[i]]
      tmp.G<-mult(XX[tmp.sis,,drop=FALSE])
      tmp.G<-tmp.G*G[tmp.anc,,drop=FALSE]
      tmp.max<-max(tmp.G)
      tmp.G<-tmp.G/tmp.max
      G[i,]<-prop.liks(tmp.G,rate,jump,sig2,elen[i])
    }
    PP<-G*X
    #divided by dx to integrate to 1
    sums<-Reduce('+',asplit(PP,2))*dx
    PP<-sweep(PP,1,sums,'/')
    attr(PP,'xpts')<-xpts
    PP
  }

  map.fun<-function(par,nt=100,
                    stochastic=FALSE,nsim=100){
    #par<-exp(par)
    rate<-par[1]
    jump<-par[2]
    sig2<-par[3]

    ##discretizing time##
    elen<-c(elen,0)
    #"ideal" time interval based on nt
    dt<-hgt/nt
    #find time res along each branch
    nts<-ceiling(elen/dt)
    #actual time interval along each branch
    dts<-elen/nts
    dts[is.nan(dts)]<-0 #for 0-length branches
    nts.seq<-lapply(nts,seq_len)

    ##new output matrix list##
    new.X<-lapply(nts+1,function(ii) matrix(nrow=ii,ncol=nx))
    for(i in tips){
      new.X[[i]][nts[i]+1,]<-X[i,,drop=FALSE]
    }
    X<-new.X

    ##dirac delta DFT function##
    DFT.dir<-dirac.DFT(nx,dx,x0=xpts[1])

    ##redefine helper functions##
    prop.liks<-function(X,rate,jump,sig2,dt,nt,nt.seq,
                        forwards=FALSE){
      #levy DFT for given t
      L<-DFT.L(dt,rate,jump)
      #normal DFT for given t
      N<-DFT.N(dt,sig2)
      #convolved DFT
      LN<-L*N
      if(forwards){
        #get DFTs of conditional likelihoods
        tmp.X<-DFT.X(X[1,])
        for(i in nt.seq+1){
          tmp.X<-tmp.X*LN
          #invert convolved DFT
          X[i,]<-IDFT.X(tmp.X)
        }
      }else{
        #get DFTs of conditional likelihoods
        tmp.X<-DFT.X(X[nt+1,])
        for(i in rev(nt.seq)){
          tmp.X<-tmp.X*LN
          #invert convolved DFT
          X[i,]<-IDFT.X(tmp.X)
        }
      }
      X
    }
    sample.fxn<-function(prob){
      unlist(lapply(seq_len(nsim),function(ii) sample(nx,1,prob=prob[,ii])),
             use.names=FALSE)
    }
    stochastic.prop.fw<-function(x0,X,rate,jump,sig2,dt,nt,nt.seq){
      out<-matrix(nrow=nt+1,ncol=nsim)
      out[1,]<-x0
      tmp.P<-matrix(nrow=nx,ncol=nsim)
      #levy DFT for given t
      L<-DFT.L(dt,rate,jump)
      #normal DFT for given t
      N<-DFT.N(dt,sig2)
      #convolved DFT
      LN<-L*N
      for(i in nt.seq+1){
        #get PDF by inverting convolved DFTs and multiplying with "backwards" conditional liks
        for(j in seq_len(nsim)){
          tmp.P[,j]<-IDFT.X(DFT.dir(x0[j])*LN)*X[i,,drop=FALSE]
        }
        #eliminate negative probs
        tmp.P[tmp.P<0]<-0
        #sample based on PDF
        out[i,]<-x0<-xpts[sample.fxn(tmp.P)]
      }
      out
    }
    comb.liks<-function(X){
      #mutliply conditional likelihoods to get conditional likelihoods at focal node
      X<-Reduce('*',X)
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

    ##pruning recursion##
    for(i in prune.seq){
      tmp.des<-des[[i]]
      for(j in tmp.des){
        X[[j]]<-prop.liks(X[[j]],rate,jump,sig2,dts[j],nts[j],nts.seq[[j]])
      }
      tmp<-comb.liks(lapply(X[tmp.des],function(ii) ii[1,,drop=FALSE]))
      if(is.null(tmp)){
        stop('underflow')
      }
      X[[i]][nts[i]+1,]<-tmp[[1]]
    }

    ##cladewise recursion##
    if(stochastic){
      out<-vector('list',nedge)
      tmp.P<-X[[nedge]]
      tmp.P[tmp.P<0]<-0
      out[[nedge]]<-matrix(sample(xpts,nsim,prob=tmp.P),1,nsim)
      attr(out[[nedge]],'tpts')<-0
      for(i in seq.int(nedge-1,1)){
        tmp.anc<-anc[i]
        x0<-out[[tmp.anc]][nts[tmp.anc]+1,]
        out[[i]]<-stochastic.prop.fw(x0,X[[i]],rate,jump,sig2,dts[i],nts[i],nts.seq[[i]])
        attr(out[[i]],'tpts')<-c(0,nts.seq[[i]])*dts[i]+attr(out[[tmp.anc]],'tpts')[nts[tmp.anc]+1]
      }
      out
    }else{
      PP<-G<-X
      G[[nedge]][1,]<-1
      PP[[nedge]]<-G[[nedge]]*X[[nedge]] #technically don't have to do this, but do have initialize G to have all the 1s...
      PP[[nedge]][1,]<-PP[[nedge]]/(dx*sum(PP[[nedge]]))
      attr(PP[[nedge]],'tpts')<-0
      for(i in seq.int(nedge-1,1)){
        tmp.anc<-anc[i]
        tmp.sis<-sis[[i]]
        tmp.G<-Reduce('*',
                      c(list(G[[tmp.anc]][nts[tmp.anc]+1,,drop=FALSE]),
                        lapply(X[tmp.sis],function(ii) ii[1,,drop=FALSE])))
        tmp.max<-max(tmp.G)
        tmp.G<-tmp.G/tmp.max
        G[[i]][1,]<-tmp.G
        G[[i]]<-prop.liks(G[[i]],rate,jump,sig2,dts[i],nts[i],nts.seq[[i]],forwards=TRUE)
        PP[[i]]<-G[[i]]*X[[i]]
        #divided by dx to integrate to 1
        sums<-Reduce('+',asplit(PP[[i]],2))*dx
        PP[[i]]<-sweep(PP[[i]],1,sums,'/')
        attr(PP[[i]],'tpts')<-c(0,nts.seq[[i]])*dts[i]+attr(PP[[tmp.anc]],'tpts')[nts[tmp.anc]+1]
      }
      attr(PP,'xpts')<-xpts
      PP
    }
  }

  out<-list('lik.fun'=lik.fun,'recon.fun'=recon.fun,'map.fun'=map.fun)
  out
}
