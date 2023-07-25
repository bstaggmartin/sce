#to-do:
# - figure out convenient way to "shut off" particular dynamics, like OU/jumps; both for simulating and fitting
# - do not optimize over rate bounds directly--instead pick reasonable bounds and go with it?
# - figure out normalization
### - likelihood increases when you try to normalize with dr, and decreases (less severely though) when you don't
### - I have to think about this more...maybe something to do with root integration?
### - Doesn't seem like it
### - Not too promising yet, fixed some mistakes, but log-likelihood increases substantially with greater nr
### - I think I figured it out, but there's some serious issues still...
### - Need to deal more formally with marginalizing over rates
### - Doesn't jump around as much, but log-lik decreases by 5-15 units every time nr is doubled...
### - Must be some weird over/underflow issue?
### - Oh boy, it does seem to "converge" in the right way, but it took like 256 nr!
### - How do I deal with that? The correspond state space is quite large...
### - Would it be easier to IFFT and FFT at that rate? Quite possibly, unfortunately...
### - Yeah, once you reach sufficient rate resolution, log-liks stabilize! So you just need to figure out a better means of integrating over this gigantic state space...
### - nt has very little effect once the log-liks stab7ilize
### - in the 75-tip example, nx had to be >=256 and nr >=64
### - Of course, hard to say how generalizable that is...
# - do I need to adjust dt for each edge?
# - increase lik by only ~1 log-lik unit using 75 tip data
#might be worth it?
#could also define some tolerance for the above issue

get.Q.fxn<-function(nx){
  base<-matrix(0,nx,nx)
  rows<-.row(c(nx,nx))
  cols<-.col(c(nx,nx))
  diag.inds<-rows==cols
  lo.off.diag.inds<-cols==rows-1
  up.off.diag.inds<-rows==cols-1
  #setting up stuff for getting approximate jump process CTMC Q matrix
  len<-2*nx+1
  dirc.base<-matrix(-1i*seq(0,2*pi,length.out=len)[-len],
                    nx,2*nx,byrow=TRUE)
  dirc.base<-sweep(dirc.base,1,seq_len(nx)-1,'*')
  norm.base<- -0.5*seq(0,pi,length.out=nx+1)^2
  norm.base<-c(norm.base,rev(norm.base)[-c(1,nx+1)])
  inds<-seq_len(nx)
  rev.inds<-rev(inds+nx)
  plan<-planFFT(2*nx)
  get.jump.Q<-function(dx,lambda,gam2){
    out<-exp(sweep(dirc.base,2,gam2/dx^2*norm.base,'+'))
    out<-Re(t(apply(out,1,IFFT,plan=plan)))
    out<-out[,inds,drop=FALSE]+out[,rev.inds,drop=FALSE]
    out[diag.inds]<-0
    out<-sweep(out,1,lambda/.rowSums(out,nx,nx),'*')
    out
  }
  xpts<-seq_len(nx)-1
  #maybe get rid of mu in favor of OU? Doesn't really make sense to have both...
  #I think this is correct, but it might be good idea to verify!
  OU.mult<-function(x0,dx,theta,alpha,sig2){
    tmp<-dx*alpha/sig2*(xpts*dx+x0+dx/2-theta)
    list(exp(tmp[-nx]),exp(-tmp[-1]))
  }
  function(x0,dx,theta,alpha,sig2,lambda,gam2){
    tmp<-OU.mult(x0,dx,theta,alpha,sig2)
    base[lo.off.diag.inds]<-tmp[[1]]*sig2/(2*dx^2)#+min(mu/dx,0)
    base[up.off.diag.inds]<-tmp[[2]]*sig2/(2*dx^2)#+max(mu/dx,0)
    Q<-base+get.jump.Q(dx,lambda,gam2)
    Q[diag.inds]<- -.rowSums(Q,nx,nx)
    Q
  }
}

#'@export
make.evorates<-function(tree,
                        exfac=0.5,nx=1024,nr=32,nt=100,
                        tip.sig=0,
                        integrate.root=TRUE){
  ##tree topology info##
  list2env(sce:::.get.tree.topo.info(tree,pruning=TRUE),envir=environment())
  #reorder continuous character
  cont<-cont[tips.ord]

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

  ##discretizing trait##
  #make a DFT plan for trait
  plan<-planFFT(2*nx)
  #find min and max considered values for trait
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
  #function for constructing DFTs for trait evolution
  base<-matrix(-0.5*sce:::.get.base(nx,dx),2*nx,nr)
  DFT.norms<-function(t,sig2s){
    exp(sweep(base,2,t*sig2s,'*'))
  }

  ##prepping rates##
  get.Q<-get.Q.fxn(nr)

  ##formatting trait data into joint array##
  #make array --> dim1 = trait value, dim2 = edges
  #each cell represents partial conditional likelihood of observed data given that the...
  #...descendant node of a particular edge exhibits a particular rate and trait
  X<-matrix(0,nrow=nx,ncol=nedge)
  scalar.init<- -length(tips)*log(dx) #think you need nr too?
  if(!tip.sig){
    inds<-findInterval(cont,xpts)
    inds<-inds+round((cont-xpts[inds])/dx)
    inds<-cbind(inds,tips)
    X[inds]<-1
  }else{
    X[,tips]<-dnorm(xpts,rep(cont,each=nx),tip.sig)
    X[,tips]<-sweep(X[,tips,drop=FALSE],2,.colSums(X[,tips,drop=FALSE],nx,length(tips)),'/')
    tmp.maxes<-apply(X[,tips,drop=FALSE],2,max)
    scalar.init<-scalar.init+sum(log(tmp.maxes))
    X[,tips]<-sweep(X[,tips,drop=FALSE],2,tmp.maxes,'/')
  }
  X<-array(X,c(nx,nedge,nr))
  XX<-array(0,c(nx,nedge,nr))
  # scalar.init<-scalar.init-1 #no idea where extra -1 comes from

  ##helper functions##
  parse.par<-function(par){
    par[-c(1,3)]<-exp(par[-c(1,3)])
    dr<-par[2]/nr
    get.Q(par[1],dr,par[3],par[4],par[5],par[6],par[7])
  }
  #about as fast as can be --> just adds ~20ish microseconds for pasting and type conversion
  mult<-function(XX){
    eval(str2lang(paste0('XX[,',seq_len(dim(XX)[2]),',,drop=FALSE]',collapse='*')))
  }
  prop.liks<-function(X,Q,t,sig2s,nt){
    #get DFTs of conditional likelihoods
    X<-apply(X,2,DFT.X)
    #get probability transition matrix for time interval
    P<-expm(t*Q)
    #scale normal DFTs for each state given sig2 and time interval
    N<-DFT.norms(t,sig2s)
    for(t in nt){
      X<-X*N
      X<-X%*%P
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
      lik<-.rowSums(lik,nx,nr)
      log(sum(lik^2))-log(sum(lik))+scalar#-log(nr)
    }
  }else{
    root.calc<-function(lik,scalar){
      lik[lik<0]<-0
      lik<-.rowSums(lik,nx,nr)
      log(max(lik))+scalar#-log(nr)
    }
  }

  #actual function to output
  lik.fun<-function(par){
    Q<-parse.par(par)
    sig2s<-exp(seq(par[1],par[1]+exp(par[2]),length.out=nr))
    #scalar --> keeps track of how cond likelihoods are scaled to prevent underflow
    dr<-log(sig2s[2]/sig2s[1])
    scalar<-scalar.init
    # P<-expm(dt*Q)
    # N<-DFT.norms(dt,sig2s)
    for(i in prune.seq){
      tmp.des<-des[[i]]
      for(j in tmp.des){
        XX[,j,]<-prop.liks(X[,j,],Q,dts[j],sig2s,nts[[j]])
      }
      tmp<-comb.liks(XX[,tmp.des,,drop=FALSE])
      if(is.null(tmp)){
        return(Inf)
      }
      X[,i,]<-tmp[[1]]
      scalar<-scalar+tmp[[2]]
    }
    -root.calc(X[,i,],scalar)
  }

  lik.fun(c(-10,log(20),-4,log(1),log(1),log(10),log(2))) #seems quite difficult to optimize over max rate bounds...
  #but then what to pick by default?
  #theta plus/minus 10*sig2?
  #other unsolved issues--very sensitive to nr
  #think there's something I'm missing in terms of normalizing constants
  #but function does seem to "peak" in the right way, which is exciting!

  out<-list('lik.fun'=lik.fun)
  out
}
