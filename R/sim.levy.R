.conv2simmap<-function(x){
  maps<-lapply(as.list(x[["edge.length"]]),setNames,nm="state_1")
  node.states<-x[["edge"]]
  mapped.edge<-matrix(x[["edge.length"]],Nedge(x),1,
                      dimnames=list(paste(node.states[,1],node.states[,2],sep=","),state="state_1"))
  node.states[]<-"state_1"
  states<-setNames(rep("1",Ntip(x)),x[["tip.label"]])
  x[c("maps","mapped.edge","node.states","states")]<-list(maps,mapped.edge,node.states,states)
  class(x)<-c("simmap",class(x))
  x
}

sim.sce<-function(tree,
                  nsims=NULL,
                  nstates=NULL,
                  disc.mod=NULL,
                  cont.mod="JN",
                  rates=1,
                  drifts=0,
                  freqs=0,
                  jumps=0,
                  skews=0,
                  disc.root=NULL,
                  cont.root=0,
                  internal=FALSE){

  #coerce any single trees to list of trees
  if(inherits(tree,"phylo")){
    tree<-list(tree)
    class(tree)<-"multiPhylo"
  }
  if(!inherits(tree,"multiPhylo")){
    stop("tree must be either a phylo, simmap, multiPhylo, or multiSimmap object")
  }

  #determining states
  #disc.mod can either be NULL or a numeric matrix/array
  #If NULL, defaults to ER if trees aren't simmaps, but just uses simmaps otherwise
  #Set nstates accordingly
  simmaps<-unlist(lapply(tree,inherits,what="simmap"))
  if(is.null(Q)&any(!simmaps)){
    Q<-1
  }
  if(!is.null(Q)&!is.numeric(Q)){
    stop("Q must either be NULL or a numeric vector/matrix/array")
  }
  if(is.null(nstates)){
    if(is.null(Q)){
      states<-sort(unique(unlist(lapply(tree,function(ii) unlist(lapply(ii[["maps"]],names),use.names=FALSE)))))
      nstates<-length(states)
    }else{
      tmp<-dim(Q)
      if(is.null(tmp)){
        states<-names(Q)

      }else{
        states<-dimnames(Q)[1:2]

      }
      states<-sort(unique(c(names(Q),unlist(dimnames(Q),use.names=FALSE))))

    }
  }
  #After getting states and nstates, format Q matrix properly
  #After formatting Q matrix, coerce any non-simmaps to simmaps




  states<-sort(unique(unlist(lapply(tree,function(ii) unlist(lapply(ii[["maps"]],names),use.names=FALSE)))))
  if(any(!simmaps)){
    if(is.numeric(Q)){
      #convert to array
      dims<-dim(Q)
      if(is.null(dims)) dims<-length(Q)
      dims<-unlist(lapply(seq_len(3),function(ii) if(is.na(dims[ii])) 1 else dims[ii]))
      dimnms<-dimnames(Q)
      if(is.null(dimnms)) dimnms<-list(names(Q))
      dimnms<-lapply(seq_len(2),function(ii) if(is.null(dimnms[ii][[1]])) rep(NA,dims[ii]) else dimnms[[ii]])
      if(is.null(dimnms[3][[1]])) dimnms[3]<-list(NULL)
      Q<-array(as.vector(Q),dims,dimnms)
    }
  }else if(!is.null(Q)){
    warning("Basing discrete regimes off of simmaps")
    Q<-NULL
  }



  if(any(!simmaps)&is.numeric(Q)){


  }else{
    Q<-NULL
  }


  if(is.null(Q)){
    warning("Q matrix unprovided and defaulted to equal rates of 1")
    Q<-array(1,c(nstates,nstates,nsims),c(rep(list(seq_len(nstates)),2),list(NULL)))
    for(i in seq_len(nstates)){
      Q[i,i,]<-1-nstates
    }
  }else if(is.numeric(Q)){
    #convert to array
    dims<-dim(Q)
    if(is.null(dims)) dims<-length(Q)
    dims<-unlist(lapply(seq_len(3),function(ii) if(is.na(dims[ii])) 1 else dims[ii]))
    dimnms<-dimnames(Q)
    if(is.null(dimnms)) dimnms<-list(names(Q))
    dimnms<-lapply(seq_len(2),function(ii) if(is.null(dimnms[ii][[1]])) rep(NA,dims[ii]) else dimnms[[ii]])
    if(is.null(dimnms[3][[1]])) dimnms[3]<-list(NULL)
    #naming is kinda tricky...


    nas<-lapply(dimnms[1:2],is.na)



    dimnms[[1]][nas[[1]]&!nas[[2]]]<-dimnms[[2]][nas[[1]]&!nas[[2]]]
    dimnms[[2]][!nas[[1]]&nas[[2]]]<-dimnms[[1]][!nas[[1]]&nas[[2]]]
    Q<-array(as.vector(Q),dims,dimnms)


    Q<-Q[seq_len(dims[1])[seq_len(nstates)],
         seq_len(dims[2])[seq_len(nstates)],
         seq_len(dims[3])[seq_len(nsims)],
         drop=FALSE]
    dimnms<-dimnames(Q)
    nas<-lapply(dimnms[1:2],is.na)
    dimnms[[1]][nas[[1]]&!nas[[2]]]<-dimnms[[2]][nas[[1]]&!nas[[2]]]
    dimnms[[2]][!nas[[1]]&nas[[2]]]<-dimnms[[1]][!nas[[1]]&nas[[2]]]
    nas<-lapply(dimnms[1:2],is.na)
    dimnms[[1]][nas[[1]]]<-


  }else{
    stop("Q must be a numeric vector, matrix, or array")
  }



  #formatting trees...
  for(i in seq_along(tree)){
    if(!inherits(tree[[i]],"simmap")){
      if(nstates==1){
        tree[[i]]<-.conv2simmap(tree[[i]])
      }else{
        #simulate based on Q matrix
      }
    }
  }




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


