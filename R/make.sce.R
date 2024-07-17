# 2/20 observation: at least in the 2-state case, the difference in drifts between states seem identifiable
# but the overall drift is not!
# kinda what you expected...
# there might be some kind of exception based on the distribution of states though...
# if you pursue this further, a simple solution might be to fix the drift of state 1 to 0
# not sure about skews...
# note, however, that the difference does seem to often be overestimated...
# at least when there is no drift difference
# seems to just be a consequence of the dataset--just hard to estimate difference in drifts
# the drift unidentifiability seems to extend to the 3 state case!
# maybe diffs in rate could make it more identifiable...
# not really...
# all this is to say: there's something here, but it's extremely difficult to identify drift params

# some initial testing: diversitree slice sampling is slow, but works well
# ML is nice and fast
# likelihood and param estimates stay similar even as res and bounds change
# bounds can be more like 0.1 (at least for BM) with no consequences
# subplex via nloptr seems to work quite well for optimization

# library(ape)
# library(fftw)
# set.seed(123)
# tree<-phytools::pbtree(n=100,scale=1)
# test<-sim.sce(tree,
#               nstates=2,
#               disc.mods=c(1,1),
#               cont.mods="JN",
#               rates=c(1,9),
#               drifts=c(0,0),
#               freqs=c(0,0),
#               jumps=c(0,0),
#               skews=c(0,0),
#               disc.roots=1)
# phytools::phenogram(get.params(test,"tree")[[1]][[1]],test$cont[,1],ftype="off")
# cont<-test$cont[,1]
# disc<-test$disc[,1]
# disc.mods=c("ARD","SYM","ER");cont.mods="JN";rates="SD";drifts=0;freqs=0
# jumps=0;skews=0;res=1024;bds.exp.fac=0.01;integrate.root=TRUE;tip.sig=0

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
make.sce<-function(tree,disc,cont,
                   cont.se=diff(range(cont,na.rm=TRUE))/100,
                   disc.mods=c("ARD","SYM","ER"),
                   cont.mods="JN",
                   rates=c("SD","EQ"),
                   drifts=0,
                   freqs=0,
                   jumps=0,
                   skews=0,
                   res=1024,bds.exp.fac=0.5,
                   integrate.root=TRUE,
                   cores=getOption("mc.cores",1L)){
  ##tree topology info##
  list2env(.get.tree.topo.info(tree,pruning=TRUE),envir=environment())
  #intial continuous data processing
  cont<-cont[tips.ord]
  #initial discrete data processing
  if(is.matrix(disc)){
    disc<-disc[tips.ord,]
    #need to add more checks
    key<-sort(colnames(disc))
    k<-length(key)
    disc<-disc[,match(key,colnames(disc)),drop=FALSE]
    disc<-t(disc)
    nas<-is.na(disc)
    n.nas<-colSums(nas)
    inds<-n.nas>0
    n.nas<-n.nas[inds]
    tmp.sums<-colSums(disc[,inds,drop=FALSE],na.rm=TRUE)
    disc[nas]<-rep(max(0,1-tmp.sums)/n.nas,n.nas)
    disc<-t(disc)
  }else{
    disc<-disc[tips.ord]
    nas<-is.na(disc)
    disc<-strsplit(disc,'&')
    key<-sort(unique(unlist(disc,use.names=FALSE)))
    k<-length(key)
    disc<-lapply(disc,function(ii) match(ii,key))
    disc[nas]<-list(seq_len(k))
    disc<-do.call(rbind,lapply(disc,tabulate,nbin=k))
  }
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
    Q.template<-disc.mods

    #making sure rows/columns of custom Q.template are named, while respecting any given names...
    # ...and coercing Q.template into a square matrix in the process.
    #Also ensures parameter indices make sense and defaults unspecified parameters to 0
    dimnms<-dimnames(Q.template)
    rownms<-rownames(Q.template)
    if(is.null(rownms)) rownms<-rep(NA,nrow(Q.template))
    rownms[!nzchar(rownms)]<-NA
    colnms<-colnames(Q.template)
    if(is.null(colnms)) colnms<-rep(NA,ncol(Q.template))
    colnms[!nzchar(colnms)]<-NA
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



  ##discretizing continuous trait##
  #find min and max considered values for continuous trait
  xlim<-range(cont,na.rm=TRUE)+
    c(-1,1)*bds.exp.fac*diff(range(cont,na.rm=TRUE))
  #get the breaks and find the inter-break distance
  xpts<-seq(xlim[1],xlim[2],length.out=res)
  dx<-diff(xpts[c(1,2)])

  #padding, backwards/forwards fft, etc., all handled by C++ now

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

  #need better warning messages...
  if(is.null(names(cont.se))){
    names(cont.se)<-rep("",length(cont.se))
  }
  probs<-which(!(names(cont.se)%in%tree[["tip.label"]]))
  matches<-match(tree[["tip.label"]],names(cont.se))
  nas<-is.na(matches)
  if(sum(nas)>0&length(probs)>0){
    matches[nas]<-rep(probs,length.out=sum(nas))
  }
  cont.se<-setNames(cont.se[matches],tree[["tip.label"]])[tips.ord]

  if(any(is.na(cont.se))){
    est.se<-TRUE
    unfixed.tips<-tips[!is.na(cont)&is.na(cont.se)]
  }else{
    est.se<-FALSE
  }

  #seems like cont.se should be at least 2*dx after all...
  #3*dx to be safe...
  #Lead to unexpected changes in likelihood as resolution changed!
  #Instead best to just returning warnings about it
  low.sigs<-(dx-cont.se)>1e-15
  if(any(low.sigs,na.rm=TRUE)){
    warning("Standard errors for some tips are very small; to ensure numerical stability, consider increasing res or standard errors. Ideally, standard errors should be no smaller than the grid resolution, which is given by (1+2*exp.bds.fac)/(res-1)*diff(range(cont))")
  }
  cont.var<-cont.se^2

  #figured out the padding stuff --> x0 should be interval midpoint
  init.dists<-get.DFTs(res,dx,c("NO","DI"),x0=(xpts[1]+xpts[res])/2)
  unfixed.init.dists<-unfixed.NO.DFT(res,dx,x0=(xpts[1]+xpts[res])/2)

  #could definitely precompute some stuff to speed up here...
  for(i in seq_along(tips)){
    X[tips[i],,]<-
      if(is.na(cont[i])){
        cbind(disc[i,]+0i,
              matrix(vector("complex",k*(2*res-1)),k,2*res-1))
      }else if(is.na(cont.var[i])){
        rep(unfixed.init.dists[["base"]](cont[i]),each=k)+log(disc[i,])
      }else if(cont.var[i]){
        rep(init.dists[[1]](cont[i],cont.var[i]),each=k)*disc[i,]
      }else{
        rep(init.dists[[2]](cont[i]),each=k)*disc[i,] #technically unnecessary now, but oh well...
      }
  }

  ##prepping inputs for C++ fun##
  des_n<-lengths(des)
  des_pos<-c(0,cumsum(des_n)[-length(des)])
  c_des<-unlist(des,use.names=FALSE)-1
  c_prune<-prune.seq-1
  c_tips<-tips-1

  ##setting up parameter parsing functionality##
  conv.dists<-get.DFTs(res,dx,c("JN","VG","NIG","BM"),
                       char.exp=TRUE)
  par.ind<-all.params!=0
  cols<-.col(dim(all.params))
  par.pos<-(cols<(k+2)|cols==(k+3)|cols==(k+4))&par.ind
  par.pol<-(cols==(k+2)|cols==(k+5))&par.ind
  par.seq<-all.params[par.ind]
  par.q.seq<-seq_len(max(all.params[,seq_len(k)]))

  lik<-function(par){
    if(est.se){
      X[unfixed.tips,,]<-
        exp(sweep(X[unfixed.tips,,,drop=FALSE],
                  3,
                  unfixed.init.dists[["modder"]](exp(2*par[length(par)])),
                  "+",check.margin=FALSE))
      par<-par[-length(par)]
    }

    all.params[par.ind]<-par[par.seq]
    all.params[par.pos]<-exp(all.params[par.pos])
    all.params[par.pol]<- -all.params[par.pol]
    tmp.Q<-all.params[,seq_len(k),drop=FALSE]
    R<-array(tmp.Q,c(k,k,2*res))
    R[dg]<- -.rowSums(tmp.Q,k,k)+
      do.call(rbind,lapply(seq_len(k),
                           function(ii)
                             conv.dists[[cont.mods[ii]]](all.params[ii,k+3],
                                                         all.params[ii,k+4],
                                                         all.params[ii,k+5])+
                             conv.dists[[4]](all.params[ii,k+1],
                                             all.params[ii,k+2])))

    -sce_lik(R,X,c_des,des_pos,des_n,c_prune,elen,scalar.init)
  }

  #parallel test
  if(cores>1){
    cores<-min(cores,max(all.params)+if(est.se) 1 else 0)
  }
  if(cores>1){
    cl<-makeCluster(cores)
    clusterExport(cl,
                  varlist=c("lik",
                            "X",
                            "est.se",
                            if(est.se) "unfixed.tips",
                            "init.dists",
                            if(est.se) "unfixed.init.dists",
                            "all.params",
                            "par.ind",
                            "par.seq",
                            "par.pos",
                            "par.pol",
                            "k",
                            "res",
                            "conv.dists",
                            "cont.mods",
                            "c_des",
                            "des_pos",
                            "des_n",
                            "c_prune",
                            "elen",
                            "scalar.init"),
                  envir=environment())
    clusterEvalQ(cl,{library(sce)})
    grad_lik<-function(par,step=1e-4){
      obj<-lik(par)
      hh<-step*abs(par)
      hh[abs(par)<1]<-step
      grad<-holder<-numeric(length(hh))
      pars<-lapply(seq_along(hh),function(ii) {par[ii]<-par[ii]+hh[ii];par})
      list("objective"=obj,
           "gradient"=(unlist(parLapply(cl,pars,function(ii) lik(ii)),use.names=FALSE)-obj)/hh)
    }
  }else{
    grad_lik<-function(par,step=1e-4){
      obj<-lik(par)
      hh<-step*abs(par)
      hh[abs(par)<1]<-step
      grad<-holder<-numeric(length(hh))
      for(i in seq_along(hh)){
        holder[i]<-hh[i]
        grad[i]<-(lik(par+holder)-obj)/hh[i]
        holder[i]<-0
      }
      list("objective"=obj,"gradient"=grad)
    }
  }

  get.quants<-function(x,probs){
    cdf<-cumsum(x)
    inds<-findInterval(probs,cdf,all.inside=TRUE)
    (probs-cdf[inds])*dx/(cdf[inds+1]-cdf[inds])+xpts[inds]
  }

  recon<-function(par,conf.lev=0.95,probs=NULL){
    if(est.se){
      X[unfixed.tips,,]<-
        exp(sweep(X[unfixed.tips,,,drop=FALSE],
                  3,
                  unfixed.init.dists[["modder"]](exp(2*par[length(par)])),
                  "+",check.margin=FALSE))
      par<-par[-length(par)]
    }

    all.params[par.ind]<-par[par.seq]
    all.params[par.pos]<-exp(all.params[par.pos])
    tmp.Q<-all.params[,seq_len(k),drop=FALSE]
    R<-array(tmp.Q,c(k,k,2*res))
    bw_R<-aperm(R,c(2,1,3))
    bw_R[dg]<- -.rowSums(tmp.Q,k,k)+
      do.call(rbind,lapply(seq_len(k),
                           function(ii)
                             conv.dists[[cont.mods[ii]]](all.params[ii,k+3],
                                                         all.params[ii,k+4],
                                                         all.params[ii,k+5])+
                             conv.dists[[4]](all.params[ii,k+1],
                                             all.params[ii,k+2])))
    all.params[par.pol]<- -all.params[par.pol]
    R[dg]<- -.rowSums(tmp.Q,k,k)+
      do.call(rbind,lapply(seq_len(k),
                           function(ii)
                             conv.dists[[cont.mods[ii]]](all.params[ii,k+3],
                                                         all.params[ii,k+4],
                                                         all.params[ii,k+5])+
                             conv.dists[[4]](all.params[ii,k+1],
                                             all.params[ii,k+2])))
    out<-sce_rec(R,bw_R,X,c_des,des_pos,des_n,c_prune,elen,c_tips)
    tmp.nms<-as.character(nodes[,2])
    tmp.nms[tips]<-tips.ord
    tmp.ord<-order(nodes[,2])
    tmp.nms<-tmp.nms[tmp.ord]
    out<-out[tmp.ord,,,drop=FALSE]
    out[out<1e-16]<-1e-16
    out<-out/apply(out,1,sum)
    dimnames(out)<-list("node"=tmp.nms,
                        "disc"=key,
                        "cont"=format(round(xpts,ceiling(-log(dx,10)))))
    marg.probs<-apply(out,c(1,2),sum)
    dimnames(marg.probs)<-list("node"=tmp.nms,
                               "disc"=key)
    if(is.null(probs)){
      probs<-(1+c(0,-1,1)*conf.lev)/2
    }
    probs[probs<0]<-0
    probs[probs>1]<-1
    marg.quants<-t(apply(out,1,function(ii) get.quants(colSums(ii),probs)))
    dimnames(marg.quants)<-list("node"=tmp.nms,
                                "quantile"=paste0(format(round(probs*100,digits=1),nsmall=1),"%"))
    attr(out,"xpts")<-xpts
    list(disc=marg.probs,
         cont=marg.quants,
         joint=out)
  }

  # out<-list('lik.fun'=lik.fun,'recon.fun'=recon.fun)
  # attr(out,'key')<-key
  # attr(out,'Q.template')<-Q.template
  dimnames(all.params)<-list("disc"=key,
                             c(key,"rates","drifts","freqs","jumps","skews"))
  list("lik"=lik,"recon"=recon,
       "se_unfixed"=est.se,
       "param_key"=all.params,
       "cont_models"=c("JN","VG","NIG")[cont.mods],
       "grad_lik"=grad_lik)
}

#old commented up versions

# lik<-function(par){
#   #could probably be optimized, but only takes ~1 millisecond for res of 1024 and 4 states
#   #so fast enough for now
#   all.params[par.ind]<-par[par.seq]
#   all.params[par.pos]<-exp(all.params[par.pos])
#   all.params[par.pol]<- -all.params[par.pol]
#   tmp.Q<-all.params[,seq_len(k),drop=FALSE]
#   R<-array(tmp.Q,c(k,k,2*res))
#   R[dg]<- -.rowSums(tmp.Q,k,k)+
#     do.call(rbind,lapply(seq_len(k),
#                          function(ii)
#                            conv.dists[[cont.mods[ii]]](all.params[ii,k+3],
#                                                        all.params[ii,k+4],
#                                                        all.params[ii,k+5])+
#                            conv.dists[[4]](all.params[ii,k+1],
#                                            all.params[ii,k+2])))
#
#   -sce_lik(R,X,c_des,des_pos,des_n,c_prune,elen,scalar.init)
# }
#
# # get.confint<-function(x,conf.lev){
# #   cdf<-cumsum(x)
# #   lb<-0.5-conf.lev/2
# #   ub<-0.5+conf.lev/2
# #   lb.ind<-which(cdf>lb)[1]
# #   ub.ind<-which(cdf>ub)[1]
# #   c(if(lb.ind==1) xpts[1] else (lb-cdf[lb.ind-1])*dx/(cdf[lb.ind]-cdf[lb.ind-1])+xpts[lb.ind-1],
# #     (ub-cdf[ub.ind-1])*dx/(cdf[ub.ind]-cdf[ub.ind-1])+xpts[ub.ind-1])
# # }
#
# get.quants<-function(x,probs){
#   cdf<-cumsum(x)
#   inds<-findInterval(probs,cdf,all.inside=TRUE)
#   (probs-cdf[inds])*dx/(cdf[inds+1]-cdf[inds])+xpts[inds]
# }
#
# recon<-function(par,conf.lev=0.95,probs=NULL){
#   all.params[par.ind]<-par[par.seq]
#   all.params[par.pos]<-exp(all.params[par.pos])
#
#   tmp.Q<-all.params[,seq_len(k),drop=FALSE]
#   R<-array(tmp.Q,c(k,k,2*res))
#
#   bw_R<-aperm(R,c(2,1,3))
#   bw_R[dg]<- -.rowSums(tmp.Q,k,k)+
#     do.call(rbind,lapply(seq_len(k),
#                          function(ii)
#                            conv.dists[[cont.mods[ii]]](all.params[ii,k+3],
#                                                        all.params[ii,k+4],
#                                                        all.params[ii,k+5])+
#                            conv.dists[[4]](all.params[ii,k+1],
#                                            all.params[ii,k+2])))
#
#   all.params[par.pol]<- -all.params[par.pol]
#   R[dg]<- -.rowSums(tmp.Q,k,k)+
#     do.call(rbind,lapply(seq_len(k),
#                          function(ii)
#                            conv.dists[[cont.mods[ii]]](all.params[ii,k+3],
#                                                        all.params[ii,k+4],
#                                                        all.params[ii,k+5])+
#                            conv.dists[[4]](all.params[ii,k+1],
#                                            all.params[ii,k+2])))
#
#   out<-sce_rec(R,bw_R,X,c_des,des_pos,des_n,c_prune,elen,c_tips)
#
#   tmp.nms<-as.character(nodes[,2])
#   tmp.nms[tips]<-tips.ord
#   tmp.ord<-order(nodes[,2])
#   tmp.nms<-tmp.nms[tmp.ord]
#   out<-out[tmp.ord,,,drop=FALSE]
#   #standardize to have common floor?
#   out[out<1e-16]<-1e-16
#   out<-out/apply(out,1,sum)
#   dimnames(out)<-list("node"=tmp.nms,
#                       "disc"=key,
#                       "cont"=format(round(xpts,ceiling(-log(dx,10)))))
#
#   marg.probs<-apply(out,c(1,2),sum)
#   dimnames(marg.probs)<-list("node"=tmp.nms,
#                              "disc"=key)
#
#   #decided to ultimately just use custom quantile function to extract central tendencies as well
#   #means can get wonky due to bimodal distributions
#
#   if(is.null(probs)){
#     probs<-(1+c(0,-1,1)*conf.lev)/2
#   }
#   probs[probs<0]<-0
#   probs[probs>1]<-1
#   marg.quants<-t(apply(out,1,function(ii) get.quants(colSums(ii),probs)))
#   dimnames(marg.quants)<-list("node"=tmp.nms,
#                               "quantile"=paste0(format(round(probs*100,digits=1),nsmall=1),"%"))
#
#   attr(out,"xpts")<-xpts
#   list(disc=marg.probs,
#        cont=marg.quants,
#        joint=out)
#
#   #decided not to do conditional probs for now...
#   #probs make a function later to pull that info out, but right now just clutters things
#   #for not much benefit
#
#   # marg.probs<-apply(out,c(1,2),sum)
#   # colnames(marg.probs)<-paste0("Pr(",key,")")
#   # colnames(marg.probs)<-paste0("prob_",key)
#
#   #found it's best to add a dx/2 correction term I think...unsure why
#   # marg.means<-as.matrix(apply(out,1,function(ii) sum(xpts*colSums(ii))+dx/2))
#   # colnames(marg.means)<-"E(marginal)"
#   # colnames(marg.means)<-"mean"
#
#   # cond.means<-apply(out,c(1,2),function(ii) sum(xpts*ii/sum(ii))+dx/2)
#   # cond.means[is.nan(cond.means)]<-NA
#   # colnames(cond.means)<-paste0("E(",key,")")
#
#   # marg.confints<-t(apply(out,1,function(ii) get.confint(colSums(ii),conf.lev=conf.lev)))
#   # colnames(marg.confints)<-paste0(round(c(0.5-conf.lev/2,0.5+conf.lev/2),3),
#   #                                 "%(marginal)")
#   # colnames(marg.confints)<-paste0(conf.lev,"%_",c("lower","upper"))
#
#   # cond.confints<-aperm(apply(out,c(1,2),function(ii) get.confint(ii/sum(ii),conf.lev=conf.lev)),
#   #                      c(2,1,3))
#   # cond.confints<-matrix(cond.confints,dim(cond.confints)[1],prod(dim(cond.confints)[-1]))
#   # colnames(cond.confints)<-paste0(round(c(0.5-conf.lev/2,0.5+conf.lev/2),3),
#   #                                 "%(",
#   #                                 rep(key,each=2),
#   #                                 ")")
#
#   # dists.sum<-cbind(marg.probs,
#   #                  marg.means,
#   #                  marg.confints)
#   #                  cond.means,
#   #                  cond.confints)
#   # rownames(dists.sum)<-tmp.nms
#
#   # attr(out,"xpts")<-xpts
#   #
#   # list("dists_sum"=dists.sum,
#   #      "dists"=out)
# }
