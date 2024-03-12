library(ape)
library(fftw)
set.seed(123)
tree<-phytools::pbtree(n=200,scale=1)
test<-sim.sce(tree,
              nstates=4,
              disc.mods=matrix(1,4,4),
              cont.mods="JN",
              rates=c(1,0,4,0),
              drifts=c(0,0,0,0),
              freqs=c(0,3,0,3),
              jumps=c(0,3,0,9),
              skews=c(0,0,0,0),
              disc.roots="STATE_1")
phytools::phenogram(get.params(test,"tree")[[1]][[1]],test$cont[,1],ftype="off")
cont<-test$cont[,1]
disc<-test$disc[,1]
disc.mods=c("ARD","SYM","ER");cont.mods="JN";rates="ARD";drifts=0;freqs=0
jumps=0;skews=0;res=1024;bds.exp.fac=0.5;integrate.root=TRUE;tip.sig=0

#generalizing...
#disc.mod may be a Q matrix template, with diagonal entries ignored
#cont.mod is a vector of state models, defaulting to BM for every state but allowing "JN", "VG", and "NIG" as well
#par.template controls how cont.mod parameters are structured:
# eventually will look like matrices with each row corresponding to a state and each column to a cont.mod parameter
# i.e., sig2 = rate, mu = drift, lambda = freq, tau2 = jump, delta = skew
# can be coarsely controlled by passing vector of "flags": AXD for all X diff (where X is R for rate, D for drift, etc.)...
# ...or AX0 for all X set to 0 (for disabling drift, BM, etc.)
# can be more finely controlled by passing template matrix

#2/15: decided instead to do separate inputs
#now only flags are SD for state-dependent and EQ for equal
#can also pass integers to assign to specific parameters
#0 just sets parameter to 0 automatically
make.sce<-function(tree,disc,cont,cont.se=0,
                   disc.mods=c("ARD","SYM","ER"),
                   cont.mods="JN",
                   rates=c("SD","EQ"),
                   drifts=0,
                   freqs=0,
                   jumps=0,
                   skews=0,
                   res=1024,bds.exp.fac=0.5,
                   integrate.root=TRUE){
  ##tree topology info##
  list2env(.get.tree.topo.info(tree,pruning=TRUE),envir=environment())
  #intial continuous data processing
  cont<-cont[tips.ord]
  #initial discrete data processing
  disc<-disc[tips.ord]
  nas<-is.na(disc)
  disc<-strsplit(disc,'&')
  key<-sort(unique(unlist(disc,use.names=FALSE)))
  k<-length(key)
  disc<-lapply(disc,function(ii) match(ii,key))
  disc[nas]<-list(seq_len(k))
  rows<-.row(c(k,k))
  cols<-.col(c(k,k))
  lt<-rows>cols
  ut<-rows<cols
  dg<-rows==cols
  odg<-lt|ut

  #discrete model
  if(!is.matrix(disc.mods)){
    ##2/15: THIS PART SEEMS TO WORK WELL##

    if(is.numeric(disc.mods)) disc.mods<-c('ARD','SYM','ER')[disc.mods[1]]
    code<-pmatch(disc.mods[1],c('ARD','SYM','ER'),nomatch=1)
    Q.template<-matrix(0,nrow=k,ncol=k,
                       dimnames=list(key,key))
    if(code==2){
      Q.template[lt]<-seq_len((k^2-k)/2)
      Q.template[ut]<-t(Q.template)[ut]
    }else{
      Q.template[odg]<-if(code==1) seq_len(k^2-k) else 1
    }
  }else{
    ##2/15: STILL NEED TO MAKE SURE THIS ALL WORKS##

    #making sure rows/columns of custom Q.template are named, while respecting any given names...
    # ...and coercing Q.template into a square matrix in the process.
    #Also ensures parameter indices make sense and defaults unspecified parameters to 0
    dimnms<-dimnames(Q.template)
    rownms<-rownames(Q.template)
    if(is.null(rownms)) rownms<-rep(NA,nrow(Q.template))
    rownms[!nzchar]<-NA
    colnms<-colnames(Q.template)
    if(is.null(colnms)) colnms<-rep(NA,ncol(Q.template))
    colnms[!nzchar]<-NA
    colnms[is.na(colnms)]<-rownms[is.na(colnms)]
    rownms[is.na(rownms)]<-colnms[is.na(rownms)]
    rowmatches<-pmatch(rownms,key)
    inds<-!is.na(rowmatches)
    rownms[inds]<-key[rowmatches[inds]]
    colmatches<-pmatch(colnms,key)
    inds<-!is.na(colmatches)
    colnms[inds]<-key[colmatches[inds]]

    #use extraneous columns to "fill in" for any unspecified rows/columns?

    dimnames(Q.template)<-list(rownms,colnms)
    Q.template<-Q.template[key,key]
    dimnames(Q.template)<-list(key,key)
    par.vec<-sort(unique(as.vector(Q.template)))
    par.vec<-par.vec[par.vec!=0&!is.na(par.vec)]
    Q.template[]<-match(Q.template,par.vec,nomatch=0L)
  }
  # par.vec<-rep(0,k^2)
  # nz.inds<-Q.template!=0

  #continuous model
  nms<-names(cont.mods)
  if(is.null(nms)){
    nms<-rep("",length(cont.mods))
  }
  if(is.numeric(cont.mods)) cont.mods<-c("JN","VG","NIG")[cont.mods]
  code<-pmatch(cont.mods,c("JN","VG","NIG"),nomatch=1)
  prob.nms<-is.na(nms)|!nzchar(nms)|!(nms%in%key)
  cont.mods<-code[match(key,nms)]
  probs<-is.na(cont.mods)
  #need better warning messages
  if(any(probs)){
    if(any(prob.nms)){
      cont.mods[probs]<-rep(code[prob.nms],length.out=sum(probs))
    }else{
      cont.mods[probs]<-1
    }
  }


  #continuous model parameters
  #now all inputs are taken separately, but coerce into matrix for ease...
  #not super-thoroughly tested, but seems to work and be pretty intuitive...
  #currently doesn't allow for the same parameters to be applied to different things
  #(i.e., param 1 couldn't correspond to both rates and drifts)
  #fine for now; may want to elaborate on this system in the future
  cont.params<-list(rates=rates,
                    drifts=drifts,
                    freqs=freqs,
                    jumps=jumps,
                    skews=skews)
  max.par.n<-max(Q.template)
  for(i in seq_along(cont.params)){
    if(is.null(cont.params[[i]])){
      cont.params[[i]]<-setNames(rep(0,k),key)
    }else if(is.character(cont.params[[i]])){
      code<-pmatch(cont.params[[i]][1],c("SD","EQ"),nomatch=1)
      cont.params[[i]]<-setNames(if(code==1) seq_along(key) else rep(1,k),
                                 key)
    }
    if(is.numeric(cont.params[[i]])){
      nms<-names(cont.params[[i]])
      if(is.null(nms)){
        nms<-rep("",length(cont.params[[i]]))
      }
      prob.nms<-is.na(nms)|!nzchar(nms)|!(nms%in%key)
      tmp<-cont.params[[i]][match(key,nms)]
      probs<-is.na(tmp)
      #need better warning messages
      if(any(probs)){
        if(any(prob.nms)){
          tmp[probs]<-rep(cont.params[[i]][prob.nms],length.out=sum(probs))
        }else{
          tmp[probs]<-0
        }
      }
      cont.params[[i]]<-tmp
    }else{
      #make more informative
      stop("Didn't recognize continuous model parameter inputs...")
    }
    par.vec<-sort(unique(cont.params[[i]]))
    par.vec<-par.vec[par.vec!=0&!is.na(par.vec)]
    cont.params[[i]]<-match(cont.params[[i]],par.vec)+max.par.n
    cont.params[[i]][is.na(cont.params[[i]])]<-0L
    max.par.n<-max(max.par.n,cont.params[[i]])
  }
  all.params<-do.call(cbind,c(list(Q.template),cont.params))
  par.inds<-match(seq_len(max.par.n),all.params)



  ##discretizing continuous trait##
  #make a DFT plan
  plan<-planFFT(2*res)
  #find min and max considered values for continuous trait
  xlim<-range(cont)+c(-1,1)*bds.exp.fac*diff(range(cont))
  #get the breaks and find the inter-break distance
  xpts<-seq(xlim[1],xlim[2],length.out=res)
  dx<-diff(xpts[c(1,2)])
  #make a vector of 0s for padding
  pad.holder<-rep(0,2*res)
  nonpad.inds1<-c((3*res/2+2):(2*res),1:(res/2+1))
  nonpad.inds2<-rep(c(TRUE,FALSE),each=res)
  pad<-function(x){
    pad.holder[nonpad.inds1]<-x
    pad.holder
  }
  #function for getting/inverting DFTs of conditional likelihoods
  DFT<-function(x) FFT(pad(x),plan=plan)
  #running into some annoying numerical issues...
  #this function is an inexpensive means of "clean up"
  IDFT<-function(x){
    tmp<-Re(IFFT(x,plan=plan))
    runs<-rle(tmp>0)
    inds<-runs$values&(runs$lengths>1)
    runs$values[inds]<-FALSE
    runs$values[!inds]<-TRUE
    tmp[inverse.rle(runs)]<-1e-16
    tmp[tmp<1e-16]<-1e-16
    tmp<-tmp/sum(tmp)*Re(x[1])
    tmp[nonpad.inds1]

    #below is equiv and about the same speed
    # tmpy<-Re(IFFT(x,plan=plan))
    # tmpy[(tmpy[c(2*res,seq_len(2*res-1))]<0&tmpy[c(2:(2*res),1)]<0)|tmpy<1e-16]<-1e-16
    # tmpy<-tmpy/sum(tmpy)*Re(x[1])
    # tmpy[nonpad.inds2]
  }


  ##formatting trait data into joint array##
  #make array --> dim1 = edges, dim2 = discrete trait, dim3 = continuous trait
  #each cell represents partial conditional likelihood of observed data given that the...
  #...descendant node of a particular edge exhibits a particular discrete-continuous phenotype
  #NEW: try keeping it on frequency spectrum this time?
  #XX is an identically-formatted array used for holding intermediate calculations
  ##TODO: add support for more flexible intraspecific distributions
  ####...both in terms of continuous dists and discrete states
  ####(i.e., prior on being in 1 discrete state vs. other)
  ntips<-length(tips)
  XX<-X<-array(0,c(nedge,k,2*res))
  scalar.init<- -ntips*log(dx)
  #figured it out!
  # tmpy<-get.DFTs(res,dx,dist="NO",x0=xpts[1])
  # plot(dnorm(xpts,xpts[500]-dx/2,sqrt(10))~xpts,type="l")
  # lines(IDFT(tmpy[[1]](xpts[500]-dx/2,10))/dx~xpts,col="red")
  #no need to deal with consts anymore
  #or round!
  #figured out the padding stuff --> x0 should be interval midpoint
  init.dists<-get.DFTs(res,dx,c("NO","DI"),x0=(xpts[1]+xpts[res])/2)

  #come up with better matching functionality for tip error...
  #seems like tip.sig should be at least 2*dx after all...
  #3*dx to be safe...
  low.sigs<-tip.sig<3*dx
  if(any(low.sigs)){
    warning("For numerical stability, tip error MUST be at least three times the resolution")
    tip.sig[low.sigs]<-3*dx
  }
  tip.sig2<-rep(tip.sig^2,length.out=ntips)

  for(i in seq_along(tips)){
    X[tips[i],disc[[i]],]<-rep(
      if(tip.sig2[i])
        init.dists[[1]](cont[i],tip.sig2[i])
      else
        init.dists[[2]](cont[i]), #technically unnecessary now, but oh well...
      each=length(disc[[i]]))
  }

  #here's the basic idea for potentially integrating C:
  # - first do eigendecomp of array in C
  # - exponentiate eigenvalues in 1 3D array using R for better stability (check first to see if this issue persists across platforms)
  # - use 3D array as input in C-based pruning algorithm that takes care of all calculations
  # - will have to incorporate FFT and IFFT cleaning procedure into C
  # - remember to check their stability!
  # - never mind! just figured it out! you need to use .st() instead of .t() for complex matrices
  # - so I think the trickiest part will just be implementing your FFT cleaning procedure in C
  # - BUT the whole pruning alg can be ported over to C then!
  # - Probs just use R to calculate transition rate array?
  # - I mean, could do that in C too...just feels unnecessary
  # - Start with parameter parsing and transition rate array construciton in R
  # - Move over to C later if you feel it's needed

  #NOTE: to make functions properly "backwards", you MUST flip the sign of any drifts/skews!

  # disc.mods=matrix(c(-1,0.1,1,-0.1),2,2),
  # cont.mods="JN",
  # rates=c(1,0),
  # drifts=c(0,0),
  # freqs=c(0,3),
  # jumps=c(0,9),
  # skews=c(0,0),

  conv.dists<-get.DFTs(res,dx,c("JN","VG","NIG","BM"),
                       char.exp=TRUE)
  testy<-array(0.1,c(k,k,2*res))
  testy[dg]<- -0.1+do.call(rbind,lapply(seq_len(3),
                          function(ii)
                            conv.dists[[c(1,1,1,1)[ii]]](c(0,3,0,3)[ii],c(0,3,0,9)[ii],0)+
                            conv.dists[[4]](c(1,0,4,0)[ii],0)))
  eigs<-lapply(asplit(testy,3),eigen)
  eigs<-lapply(eigs,function(ii) setNames(c(ii,list((solve(ii$vectors)))),
                                          c("eigvals","eigvecs","inveigvecs")))
  foo<-function(t,X){
    tmp<-lapply(eigs,function(ii) ii$eigvecs%*%diag(exp(t*ii$eigvals))%*%ii$inveigvecs)
    matrix(
      unlist(lapply(seq_along(tmp),function(ii) tmp[[ii]]%*%X[,ii]),use.names=FALSE),
      k,2*res)
  }
  testy.out<-apply(foo(elen[1],X[1,,]),1,IDFT)*apply(foo(elen[1],X[1,,]),1,IDFT)
  plot((Re(fftw::IFFT(testy.out[2,]))))
  matplot(apply(testy.out,2,DFT),type="l") #Periodic boundary conditions might prove more problematic than anticipated...
  #something to work on later...more relevant for Levy processes

  eigs2<-vec_eig(testy)
  testy.out2<-vec_exp_mult(0.2,X[tips[1],,],eigs2[[1]],eigs2[[2]],eigs2[[3]])
  matplot(log(apply(testy.out2,1,IDFT)),type="l") #works!

  microbenchmark::microbenchmark(foo(0.2))
  microbenchmark::microbenchmark(vec_exp_mult(0.2,X[tips[1],,],eigs2[[1]],eigs2[[2]],eigs2[[3]]))
  #sweet! about a 40-fold speedup (though what you have above isn't necessarily the most efficient implementation...)
  #still, once you have the entire pruning algorithm in C, should be a breeze

  iii<-3
  #should work...definitely seems like a math error
  #there's something you're not getting...
  (eigs2[[2]][,,iii]%*%(diag(exp(eigs2[[1]][,iii]))%*%(eigs2[[3]][,,iii]%*%X[tips[1],,iii])))
  (eigs[[iii]][[2]]%*%(diag(exp(eigs[[iii]][[1]]))%*%(eigs[[iii]][[3]]%*%X[tips[1],,iii])))


  testy.out2[,iii]
  eigs2[[3]][,,iii]%*%X[tips[1],,iii]

  testy.out3<-lapply(seq_along(eigs),function(ii) eigs[[ii]]$inveigvecs%*%X[tips[1],,ii])

  matplot(log(apply(testy.out2,1,IDFT)),type="l") #hmmm...feels like a math error
  testy.out3<-test_exp(0.2,eigs2[[1]])
  matplot(t(testy.out3),type="l")
  testy.out4<-exp(0.2*do.call(cbind,lapply(eigs,"[[","eigvals")))
  matplot(t(testy.out4),type="l")

  testy.out3<-apply(testy.out3,2,sort)
  testy.out4<-apply(test.out4,2,sort)
  plot(Re(testy.out4[1,])~Re(testy.out3[1,]))
  #so they DO return the same results! what gives?

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

#helpful C functions
#R may be way more stable for eigendecomp
#But worth trying
#Just checked--numerical stability seems fine!
#Yeah, 50-fold increase in speed--totally worth it

Rcpp::cppFunction("
List vec_eig(arma::cx_cube x) {

  arma::cx_mat eigvals(x.n_rows, x.n_slices);
  arma::cx_vec eigval(x.n_rows);
  arma::cx_cube eigvecs(x.n_rows, x.n_cols, x.n_slices);
  arma::cx_cube inveigvecs = eigvecs;
  arma::cx_mat eigvec(x.n_rows, x.n_cols);
  // arma::uvec indices;

  for(int i=0; i<x.n_slices; i++) {
    eig_gen(eigval, eigvec, x.slice(i));
    // indices = arma::sort_index(eigval);
    eigvals.col(i) = eigval;
    eigvecs.slice(i) = eigvec;
    inveigvecs.slice(i) = arma::inv(eigvec);
  }

  return List::create(eigvals, eigvecs, inveigvecs);
}
",
depends="RcppArmadillo")


#figured it out!!!
#it seems like, for whatever reason, the complex and imaginary components get "swapped" by Rcpp!
#or something like that...
Rcpp::cppFunction("
arma::cx_mat vec_exp_mult(double t, arma::cx_mat x, arma::cx_mat eigvals, arma::cx_cube eigvecs, arma::cx_cube inveigvecs) {

  arma::cx_mat out1(eigvals.n_rows, eigvals.n_cols);
  arma::cx_mat tmp_exp(eigvals.n_rows, eigvals.n_cols);
  arma::cx_mat out2(eigvals.n_rows, eigvals.n_cols);

  // something seems to be going wrong below...
  // even if you sort the eigenvalues and such beforehand...
  // I figured it out! .t() does the conjugate transpose for complex matrices
  // You want .st() instead!
  for(int i=0; i<eigvals.n_rows; i++) {
    out1.row(i) = arma::sum(inveigvecs.row_as_mat(i).st() % x);
  }
  out1 %= arma::exp(t * eigvals);
  for(int i=0; i<eigvals.n_rows; i++) {
    out2.row(i) = arma::sum(eigvecs.row_as_mat(i).st() % out1);
  }

  return out2;
}
",
depends="RcppArmadillo")

Rcpp::cppFunction("
arma::cx_mat test_exp(double t, arma::cx_mat eigvals) {

  return arma::exp(t * eigvals);

}
",
depends="RcppArmadillo")

#
# #Will have to pass my own exponentiated eigenvalues
# #But other than that should work
# Rcpp::cppFunction("
# arma::cx_mat vec_exp_mult(double t, arma::cx_mat x, arma::cx_mat eigvals, arma::cx_cube eigvecs, arma::cx_cube inveigvecs) {
#
#   arma::cx_mat out1(eigvals.n_rows, eigvals.n_cols);
#   arma::cx_mat tmp_exp(eigvals.n_rows, eigvals.n_cols);
#   arma::cx_mat out2(eigvals.n_rows, eigvals.n_cols);
#
#   for(int i=0; i<eigvals.n_rows; i++) {
#     out1.row(i) = arma::sum(inveigvecs.row_as_mat(i).t() % x);
#   }
#   // out1 %= arma::exp(t * eigvals);
#   tmp_exp.set_real(exp(t * real(eigvals)) % cos(t * imag(eigvals)));
#   tmp_exp.set_imag(exp(t * real(eigvals)) % sin(t * imag(eigvals)));
#   out1 %= tmp_exp;
#   for(int i=0; i<eigvals.n_rows; i++) {
#     out2.row(i) = arma::sum(eigvecs.row_as_mat(i).t() % out1);
#   }
#
#   return out2;
# }
# ",
# depends="RcppArmadillo")
