
.getstates<-function(x){
  sort(
    unique(
      unlist(
        lapply(x,function(ii) unlist(lapply(ii[["maps"]],names),use.names=FALSE)),
        use.names=FALSE)
      )
    )
}

.fixvec<-function(Q,states,k,dims){
  #form basic vector
  out<-setNames(Q[states],states)
  probs<-is.na(out)
  if(any(probs)){
    #recycle as necessary from unlabelled/unmatched inputs
    nms<-names(Q)
    if(is.null(nms)){
      prob.nms<-rep(TRUE,dims)
    }else{
      prob.nms<-is.na(nms)|!nzchar(nms)|!(nms%in%states)
    }
    if(any(prob.nms)){
      out[probs]<-rep(Q[prob.nms],length.out=sum(probs))
      probs<-is.na(probs)
    }
    #fill in anything remaining with defaults
    out[probs]<-1/k
  }
  out<-tmp<-out/sum(out)
  out<-array(rep(out,each=k),c(k,k,1),
              dimnames=list(states,states,NULL))
  diag.inds<-rep(rep(c(TRUE,FALSE),c(1,k)),length.out=k^2)
  out[diag.inds]<-0
  out[diag.inds]<- -apply(out,c(1,3),sum)
  out<-out/-sum(tmp*out[diag.inds])
  out
}

.fixarr<-function(Q,states,k,dims){
  #convert to array if necessary...
  if(length(dims)==2){
    Q<-array(Q,c(dims,1),dimnames(Q))
    dims<-c(dims,1)
  }else if(length(dims)>3){
    stop("disc.mods should have 3 dimensions at max corresponding to ingoing state, outgoing state, then simulation")
  }
  dimnms<-dimnames(Q)
  #first symmetrize names to best of ability...
  rownms<-dimnms[[1]]
  if(is.null(rownms)) rownms<-rep("",dims[1])
  prob.rownms<-is.na(rownms)|!nzchar(rownms)|!(rownms%in%states)
  colnms<-dimnms[[2]]
  if(is.null(colnms)) colnms<-rep("",dims[2])
  prob.colnms<-is.na(colnms)|!nzchar(colnms)|!(colnms%in%states)
  tmp.seq<-seq_len(min(dims[1:2]))
  tmp.inds1<-which(!prob.rownms[tmp.seq]&prob.colnms[tmp.seq])
  tmp.inds2<-which(prob.rownms[tmp.seq]&!prob.colnms[tmp.seq])
  colnms[tmp.inds1]<-rownms[tmp.inds1]
  rownms[tmp.inds2]<-colnms[tmp.inds2]
  prob.colnms[tmp.inds1]<-FALSE
  prob.rownms[tmp.inds2]<-FALSE
  #then fill in missing names...
  #NOTE: you may want to do a more careful recycling procedure where unmatched inputs can be...
  #piped back to fill in missing entries in the basic array, but there isn't really an
  #...intuitive way to do this and it seems unnecessary for now
  if(any(prob.colnms)){
    #should always have problem states in this case I think, but just in case...
    prob.states<-!(states%in%colnms)
    if(any(prob.states)){
      colnms[prob.colnms]<-rep(states[prob.states],length.out=sum(prob.colnms))
    }
  }
  if(any(prob.rownms)){
    #should always have problem states in this case I think, but just in case...
    prob.states<-!(states%in%rownms)
    if(any(prob.states)){
      rownms[prob.rownms]<-rep(states[prob.states],length.out=sum(prob.rownms))
    }
  }
  dimnms[1:2]<-list(rownms,colnms)
  dimnames(Q)<-dimnms
  #form basic array
  out<-Q[match(states,rownms),match(states,colnms),,drop=FALSE]
  dimnames(out)<-list(states,states,NULL)
  #round <0 entries up to 0
  out[out<0]<-0
  #fill in missing entries with defaults
  out[is.na(out)]<-1/(k-1)
  #fix diagonal
  diag.inds<-rep(rep(c(TRUE,FALSE),c(1,k)),length.out=k^2)
  out[diag.inds]<-0
  out[diag.inds]<- -apply(out,c(1,3),sum)
  out
}

.conv2simmap<-function(x,states.vec){
  maps<-lapply(as.list(x[["edge.length"]]),setNames,nm=states.vec)
  node.states<-x[["edge"]]
  mapped.edge<-matrix(x[["edge.length"]],Nedge(x),1,
                      dimnames=list(paste(node.states[,1],node.states[,2],sep=","),state=states.vec))
  node.states[]<-states.vec
  states<-setNames(rep("1",Ntip(x)),x[["tip.label"]])
  x[c("maps","mapped.edge","node.states","states")]<-list(maps,mapped.edge,node.states,states)
  class(x)<-c("simmap",class(x))
  x
}

#for continuous model parameters
#for conveience, recycles any unmatched inputs to match the length of prob.states
#(rather than only pair up in 1-to-1 manner)
#(such logic does not seem to apply as readily to Q coercion, so I think it's fine to not do it there)
.fixmat<-function(mat,states,nstates,nice.nm){
  if(is.null(mat)) mat<-as.numeric(NA)
  if(!is.numeric(mat)){
    if(nice.nm=="cont.mods"){
      if(is.character(mat)){
        mat[]<-pmatch(mat,c("JN","VG","NIG"))
        mode(mat)<-"numeric"
      }else{
        stop(paste0(nice.nm," should be either a character vector/matrix or a numeric one (1=JN, 2=VG, 3=NIG)"))
      }
    }else{
      stop(paste0(nice.nm," should be a numeric vector/matrix"))
    }
  }
  if(length(dim(mat))<1){
    mat<-matrix(mat,length(mat),1,
                dimnames=list(names(mat),NULL))
  }
  if(length(dim(mat))>2){
    stop(paste0(nice.nm," should only have 2 dimension at max corresponding to states than simulations"))
  }
  nms<-rownames(mat)
  if(is.null(nms)) nms<-rep("",dim(mat)[1])
  prob.nms<-is.na(nms)|!nzchar(nms)|!(nms%in%states)
  if(any(prob.nms)){
    prob.states<-!(states%in%nms)
    if(any(prob.states)){
      nms[prob.nms]<-rep(states[prob.states],length.out=sum(prob.nms))
      nrem<-sum(prob.states)-sum(prob.nms)
      if(nrem>0){
        mat<-rbind(mat,mat[rep(which(prob.nms),length.out=nrem),,drop=FALSE])
        nprobs<-sum(prob.states)
        nms<-c(nms,states[prob.states][(nprobs-nrem+1):nprobs])
      }
    }
  }
  rownames(mat)<-nms
  mat<-mat[match(states,nms),,drop=FALSE]
  rownames(mat)<-states
  mat[is.na(mat)]<-if(nice.nm=="cont.mods"|nice.nm=="rates") 1 else 0
  mat
}

.gettreeinfo<-function(x,states,internal){
  edge<-x[["edge"]]
  anc<-do.call(match,asplit(edge,2))
  maps<-x[["mapped.edge"]]
  maps<-maps[,match(states,colnames(maps)),drop=FALSE]
  dimnames(maps)<-NULL
  maps[is.na(maps)]<-0
  attr(x,"order")<-NULL
  seq<-reorder(x,index.only=TRUE)
  nn<-Ntip(x)
  node.ord<-order(c(nn+1,edge[,2]))
  lens<-lengths(x[["maps"]])
  disc.root<-names(x[["maps"]][is.na(anc)][[1]])[1]
  disc<-c(disc.root,names(unlist(x[["maps"]]))[cumsum(lens)])[node.ord]
  #just ignoring node labels for now for ease...
  nms<-x[["tip.label"]]
  if(internal){
    nms<-c(nms,(nn+1):(nn+x[["Nnode"]]))
  }else{
    tmp.seq<-seq_len(nn)
    node.ord<-node.ord[tmp.seq]
    disc<-disc[tmp.seq]
  }
  list(anc=anc,maps=maps,seq=seq,node.ord=node.ord,nms=nms,disc=setNames(disc,nms),
       disc.root=disc.root,disc.mods=x[["Q"]])
}

#make better version that avoids doing anything for entries with mean and/or shape equal to 0
#when shape is equal to 0, approaches a limiting cases where values are either 0 of Inf, it seems
#not really something I can do anything with...
#so probably best to just keep skew strictly between -1 and 1 for NIG processes
.rinvgauss<-function(n,mean,shape){
  nu<-rnorm(n)
  nu.sq<-nu^2
  mean.sq<-mean^2
  meannu.sq<-mean.sq*nu.sq
  shape2<-shape*2
  x<-mean+meannu.sq/shape2-mean/shape2*sqrt(4*mean*shape*nu.sq+meannu.sq*nu.sq)
  tmp.inds<-!is.nan(x)
  nn<-sum(tmp.inds)
  inds<-vector(length=n)
  tmp<-rep(mean,length.out=n)[tmp.inds]
  inds[tmp.inds]<-runif(nn)>tmp/(tmp+x[tmp.inds])
  x[inds]<-rep(mean.sq,length.out=n)[inds]/x[inds]
  x[!tmp.inds]<-rep(shape,length.out=n)[!tmp.inds]/rchisq(n-nn,1)
  x
}

.conv2phylo<-function(x){
  x[c("maps","mapped.edge","node.states","states","Q","logL")]<-NULL
  class(x)<-"phylo"
  x
}

# base<-phytools::pbtree(n=100,scale=1)
# qq<-matrix(c(-1,1,1,-1),2,2,
#            dimnames=list(letters[1:2],letters[1:2]))
# tree<-phytools::sim.history(base,qq,nsim=20)
# disc.mods<-setNames(c(4,1,3),letters[c(2,1,3)])
# disc.roots<-"c"
# tree[c(5,7,15:20)]<-list(base)
# nsims=NULL;nstates=NULL;cont.mods=c("JN");rates=1;drifts=-2;freqs=c(0,1,10)
# jumps=9;skews=0.5;cont.roots=-10;internal=TRUE

#need to clean up, but has all the features I want...
#' @export
sim.sce<-function(trees,
                  nsims=NULL,
                  nstates=NULL,
                  disc.mods=NULL,
                  cont.mods="JN",
                  rates=1,
                  drifts=0,
                  freqs=0,
                  jumps=0,
                  skews=0,
                  disc.roots=NULL,
                  cont.roots=0,
                  internal=FALSE){

  #coerce any single trees to list of trees
  if(inherits(trees,"phylo")){
    trees<-list(trees)
    class(trees)<-"multiPhylo"
  }
  if(!inherits(trees,"multiPhylo")){
    stop("trees must be either a phylo, simmap, multiPhylo, or multiSimmap object")
  }

  ##DISCRETE MODEL OVERVIEW##
  #>Use simmaps if available
  #>Otherwise, use disc.mods, which can either be NULL or a numeric vector/matrix/array of Q matrices
  #>Vector Q matrices default to specifying the stationary dist

  #>2/6 Aside:
  #>>Alright, I got stuck on the above because many Q's can yield the same stationary dist
  #>>So it's hard to find a "natural" Q to pick...
  #>>But then I realized you can just borrow from the general time reversible model
  #>>Assuming substitution params of 1, i.e.:
  # Q<-rexp(5)
  # Q<-in.Q<-Q/sum(Q)
  # Q<-matrix(rep(Q,each=length(Q)),length(Q),length(Q))
  # diag(Q)<-0
  # diag(Q)<- -rowSums(Q)
  # Q<-Q/-sum(in.Q*diag(Q))
  # expm::expm(Q*100);in.Q
  #>>Haha yep! That works. Moving on...

  #>If disc.mods is needed but not provided, defaults to flat stationary dist
  #>If not provided, nstates is selected based on the total number of unique state names
  #>Across any provided simmaps or Q arrays
  #>Treating unlabelled entries for a given dimension as potentially unique
  #>Though unlabelled entries will default to unfound states in simmaps first...

  #checking for simmaps and basic disc.mods formatting requirements
  simmaps<-unlist(lapply(trees,inherits,what="simmap"))
  Q<-disc.mods
  if(is.null(Q)&any(!simmaps)){
    Q<-1
  }
  if(!is.null(Q)&!is.numeric(Q)){
    stop("disc.mods must either be NULL or a numeric vector/matrix/3D array")
  }
  if(!is.null(disc.roots)&!is.numeric(disc.roots)&!is.character(disc.roots)){
    stop("disc.roots must either be NULL, a numeric vector/matrix, or a character vector")
  }

  #getting states and nstates
  #probably should have better warning messages
  #think this should work, but not thoroughly tested
  if(any(simmaps)){
    simmap.states<-.getstates(trees[simmaps])
  }else{
    simmap.states<-NULL
  }
  if(!is.null(Q)){
    dims<-dim(Q)
    if(is.null(dims)){
      ndims<-1L
      nms<-list(names(Q))
      dims<-length(Q)
    }else{
      ndims<-2L
      nms<-dimnames(Q)[1:2]
    }
    nprobs<-vector("numeric",ndims)
    for(i in 1:ndims){
      if(length(nms)<i){
        nms[[i]]<-rep("",dims[i])
      }
      nms[[i]][is.na(nms[[i]])]<-""
      nprobs[i]<-sum(!nzchar(nms[[i]]))
    }
    Q.states<-unique(unlist(nms,use.names=FALSE))
    Q.states<-Q.states[nzchar(Q.states)]
    extra.states<-max(nprobs)
    if(!is.null(simmap.states)){
      extra.states<-max(0,extra.states-sum(!(simmap.states%in%Q.states)))
    }
    if(extra.states>0){
      #will cause issues if anyone specifies a name like STATE_XX
      Q.states<-c(Q.states,paste0("STATE_",1:extra.states))
    }
  }else{
    Q.states<-NULL
    extra.states<-0
  }

  #also search for any names in disc.roots and add to Q.states-->I don't think this will be an issue...
  if(!is.null(disc.roots)){
    if(is.character(disc.roots)){
      tmp.states<-unique(disc.roots)
    }
    if(is.numeric(disc.roots)){
      if(is.null(dim(disc.roots))){
        tmp.states<-unique(names(disc.roots))
      }else{
        tmp.states<-unique(dimnames(disc.roots)[[1]])
      }
    }
    tmp.states<-tmp.states[nzchar(tmp.states)&!is.na(tmp.states)]
    Q.states<-unique(c(Q.states,tmp.states))
  }

  states<-unique(c(simmap.states,Q.states))
  if(!is.null(nstates)){
    if(!is.numeric(nstates)){
      stop("nstates must be NULL or a single integer")
    }
    if(length(nstates)>1){
      warning("nstates should be a single integer; only using first entry of nstates")
      nstates<-nstates[1]
    }
    if(is.na(nstates)|nstates<max(1,length(simmap.states))){
      warning("nstates must be at least the number of states specified by trees; defaulting to lowest possible nstates value")
      nstates<-max(1,length(simmap.states))
    }
    rem.states<-nstates-length(states)
    if(rem.states>0){
      warning("nstates implies more states than what was found in inputs; adding extra states")
      states<-c(states,paste0("STATE_",1:rem.states+extra.states))
    }else if(rem.states<0){
      warning("nstates implies fewer states than what was found in inputs; ignoring excess states found in disc.mods/disc.roots")
      states<-states[1:(length(states)+rem.states)]
    }
  }
  nstates<-length(states)

  #finally coerce trees to proper format
  #may at some point want to get disc.roots and/or Qs directly from simmap objects...
  disc.roots.lookup<-disc.mods.lookup<-rep(NA,length(trees))
  nQ<-nR<-nn<-0
  if(any(!simmaps)){
    if(nstates>1){

      #coerce Q to proper format
      ndims<-length(dims)
      #probably best to just handle each case separately because vectors are handled specially now...
      #need to implement informative warning messages
      if(length(dims)<2){
        Q<-.fixvec(Q,states,nstates,dims)
      }else{
        Q<-.fixarr(Q,states,nstates,dims)
      }

      #coerce disc.roots to proper format
      if(is.null(disc.roots)){
        disc.roots<-matrix(1/nstates,nstates,1,
                          dimnames=list(states,NULL))
      }else{
        if(is.numeric(disc.roots)){
          if(length(dim(disc.roots))<1){
            if(any(disc.roots==0)|all(disc.roots<=1)|max(abs(round(disc.roots)-disc.roots))>1e-10){
              #assume probabilities
              disc.roots<-matrix(disc.roots,length(disc.roots),1,
                                dimnames=list(names(disc.roots),NULL))
            }else{
              #assume indices otherwise
              disc.roots<-states[disc.roots]
            }
          }
          if(length(dim(disc.roots)==2)){
            nms<-rownames(disc.roots)[[1]]
            if(is.null(nms)) nms<-rep("",dim(disc.roots)[1])
            prob.nms<-is.na(nms)|!nzchar(nms)|!(nms%in%states)
            prob.states<-!(states%in%nms)
            if(any(prob.states)){
              nms[prob.nms]<-rep(states[prob.states],length.out=sum(prob.nms))
            }
            rownames(disc.roots)<-nms
            disc.roots<-disc.roots[match(states,nms),,drop=FALSE]
            rownames(disc.roots)<-states
          }else{
            stop("disc.roots should have 2 dimensions at max corresponding to state then simulation")
          }
        }
        if(is.character(disc.roots)){
          if(length(dim(disc.roots))>1){
            warning("coercing disc.roots to a character vector rather than array")
          }
          disc.roots<-as.vector(disc.roots)
          tmp.root<-matrix(0,nstates,length(disc.roots),
                           dimnames=list(states,NULL))
          tmp.root[cbind(match(disc.roots,states),seq_along(disc.roots))]<-1
          disc.roots<-tmp.root
        }
        disc.roots[is.na(disc.roots)|disc.roots<0]<-0
        tmp.sums<-colSums(disc.roots)
        disc.roots<-disc.roots/rep(tmp.sums,each=nstates)
        disc.roots[,tmp.sums==0]<-1/nstates
      }

      #generate simmaps...
      nQ<-dim(Q)[3]
      nR<-dim(disc.roots)[2]
      nn<-max(nQ,nR)
      prob.inds<-which(!simmaps)
      trees[!simmaps]<-lapply(seq_along(prob.inds),
                             function(ii) phytools::sim.history(trees[[prob.inds[ii]]],
                                                                Q=t(Q[,,(ii-1)%%nQ+1]),
                                                                anc=disc.roots[,(ii-1)%%nR+1],
                                                                message=FALSE))
      disc.mods.lookup[!simmaps]<-seq_len(nQ)
      disc.roots.lookup[!simmaps]<-seq_len(nR)
    }else{
      #generate single-state simmaps...
      trees[!simmaps]<-lapply(trees[!simmaps],.conv2simmap,states.vec=states)
    }
  }

  #And now, we FINALLY have a properly formatted set of simmaps and state info
  #Probably went overboard with coercing input again, but oh well
  #Let's get to simulating the continuous character stuff!!!

  #coerce continuous model parameters...
  params<-list(cont.mods=cont.mods,
               rates=rates,
               drifts=drifts,
               freqs=freqs,
               jumps=jumps,
               skews=skews)
  params<-lapply(seq_along(params),
                 function(ii) .fixmat(params[[ii]],
                                      states,
                                      nstates,
                                      names(params)[ii]))
  lens<-unlist(lapply(params,function(ii) dim(ii)[2]),use.names=FALSE)

  #coerce cont.roots real quick...
  if(is.null(cont.roots)) cont.roots<-as.numeric(NA)
  if(!is.numeric(cont.roots)){
    stop("cont.roots should be a numeric vector")
  }
  if(length(dim(cont.roots))>2){
    warning("coercing cont.roots to a vector rather than array")
    cont.roots<-as.vector(cont.roots)
  }
  cont.roots[is.na(cont.roots)]<-0

  #final params list
  params<-c(params,
            trees=list(lapply(trees,.gettreeinfo,states=states,internal=internal)),
            cont.roots=list(cont.roots))
  lens<-c(lens,length(trees),length(cont.roots))
  if(is.null(nsims)){
    nsims<-max(lens)
  }

  .getparams<-function(ii,i){
    if(ii>6){
      params[[ii]][[(i-1)%%lens[ii]+1]]
    }else{
      params[[ii]][,(i-1)%%lens[ii]+1]
    }
  }

  #could potentially be more intelligently vectorized...
  #but, eh, not worried for now
  cont.out<-vector("list",nsims)
  for(i in seq_len(nsims)){
    tmp<-lapply(seq_along(params),.getparams,i=i)
    maps<-tmp[[7]][["maps"]]
    nedges<-length(tmp[[7]][["anc"]])

    #need to figure NIG with skew--done!
    #avoiding simulations for 0 time increments actually doesn't speed anything up
    #but I like the clarity so I'm keeping it this way
    seed<-maps
    inds<-maps>0
    nn<-.colSums(inds,nedges,nstates)
    seed[inds]<-
      matrix(
        unlist(
          lapply(seq_len(nstates),
                 function(ii)
                   switch(tmp[[1]][ii],
                          rpois(nn[ii],maps[inds[,ii],ii]*tmp[[4]][ii]),
                          #3/7: fixed a mistake below! 3rd argument should be tmp[[4]][ii], NOT 1/tmp[[4]][ii]
                          rgamma(nn[ii],maps[inds[,ii],ii]*tmp[[4]][ii],tmp[[4]][ii]),
                          #3/7: fixed
                          .rinvgauss(nn[ii],
                                     tmp[[5]][ii]*maps[inds[,ii],ii]/sqrt(1-tmp[[6]][ii]^2),
                                     (tmp[[4]][ii]*tmp[[5]][ii]*maps[inds[,ii],ii])^2))),
          use.names=FALSE
        )
      )

    #you should be able to do this for vg, but I'm not sure about NIG...
    #can do this with NIG, but I'll have to modify tmp[[6]] and tmp[[5]] specially for NIG states
    #tmp[[6]] should be tau2*delta/lambda
    #tmp[[5]] should be tau2
    #this part could potentially be precomputed
    #since reparam'ing NIG, only need tmp.skew
    #jumps are identical across processes
    tmp.jump<-tmp[[5]]
    nig<-tmp[[1]]==3
    if(any(nig)){
      tmp.jump[nig]<-1
    }
    #I belive the math here all checks out, but be vigilant about it!
    seed<-rnorm(nedges,
                seed%*%tmp[[6]]+maps%*%tmp[[3]],
                sqrt(seed%*%tmp.jump+maps%*%tmp[[2]]))

    for(j in tmp[[7]][["seq"]]){
      a<-tmp[[7]][["anc"]][j]
      if(is.na(a)){
        seed[j]<-seed[j]+tmp[[8]]
      }else{
        seed[j]<-seed[j]+seed[a]
      }
    }
    cont.out[[i]]<-setNames(c(tmp[[8]],seed)[tmp[[7]][["node.ord"]]],
                            tmp[[7]][["nms"]])
  }
  # phytools::phenogram(tree[[16]],cont.out[[16]],ftype="off")
  #should probably check for tip/node label consistency across outputs and bind into matrices if possible...
  #or even just bind into matrices with missing entries?
  #nah, I think I like only binding into matrix assuming equal length and names...
  nms<-lapply(params[[7]],"[[","nms")
  disc.out<-lapply(params[[7]],function(ii) ii[["disc"]])
  if(all(lengths(nms[-1])==length(nms[[1]]))){
    matches<-do.call(cbind,c(list(seq_along(nms[[1]])),lapply(nms[-1],match,x=nms[[1]])))
    nn<-ncol(matches)
    if(!any(is.na(matches))){
      cont.out<-do.call(cbind,lapply(seq_len(nsims),function(ii) cont.out[[ii]][matches[,(ii-1)%%nn+1]]))
      disc.out<-do.call(cbind,lapply(seq_along(disc.out),function(ii) disc.out[[ii]][matches[,ii]]))
    }
  }
  if(is.list(disc.out)){
    disc.out<-rep(disc.out,length.out=nsims)
  }else{
    disc.out<-disc.out[,rep(seq_len(ncol(disc.out)),length.out=nsims),drop=FALSE]
  }

  #make cont.mods more intuitive
  params[[1]][]<-c("JN","VG","NIG")[params[[1]]]
  #process.disc.mods stuff...
  if(nstates==1){
    disc.roots<-matrix(1,nstates,1,
                      dimnames=list(states,NULL))
    disc.roots.lookup<-rep(1,length(trees))
    Q<-array(dim=c(nstates,nstates,0),
             dimnames=list(states,states,NULL))
    mode(Q)<-"numeric"
  }else{
    nn<-sum(simmaps)
    tmp.seq<-seq_len(nn)
    extra.disc.roots<-matrix(0,nstates,nn,
                            dimnames=list(states,NULL))
    tmp.roots<-unlist(lapply(params[[7]][simmaps],function(ii) ii[["disc.root"]]))
    extra.disc.roots[cbind(match(tmp.roots,states),tmp.seq)]<-1
    disc.roots.lookup[simmaps]<-nR+tmp.seq
    disc.roots<-cbind(disc.roots,extra.disc.roots)
    #should never replace specified Q matrices in disc.mods because sim.history doesn't return a Q matrix...
    tmp.Qs<-lapply(params[[7]][simmaps],function(ii) ii[["disc.mods"]])
    avails<-lengths(tmp.Qs)>0
    nnn<-sum(avails)
    Q<-array(if(nQ+nnn>0) c(Q,unlist(tmp.Qs,use.names=FALSE)) else 1,
             c(nstates,nstates,nQ+nnn),
             list(states,states,NULL))
    disc.mods.lookup[avails]<-nQ+seq_len(nnn)
  }

  names(params)<-c("cont.mods","rates","drifts","freqs","jumps","skews","trees","cont.roots")
  #strip simmapped histories in internal is FALSE? Good idea?
  #Nah, I feel it's more like a parameter than actual data...
  # if(!internal){
  #   trees<-lapply(trees,.conv2phylo)
  #   class(trees)<-"multiPhylo"
  # }
  params[["trees"]]<-trees
  params<-c(params,disc.mods=list(Q),disc.roots=list(disc.roots))
  lookups<-do.call(cbind,c(lapply(lens,function(ii) rep(seq_len(ii),length.out=nsims)),
                           list(rep(disc.mods.lookup,length.out=nsims)),
                           list(rep(disc.roots.lookup,length.out=nsims))))
  colnames(lookups)<-names(params)
  params<-params[c(7,9,1:6,10,8)]
  lookups<-lookups[,c(7,9,1:6,10,8),drop=FALSE]

  out<-list(disc=disc.out,cont=cont.out)
  attr(out,"param_info")<-list(params=params,lookups=lookups)
  class(out)<-"sce_sim"
  out
}

# tmp.tree<-tree[[5]]
# tmpy<-seed[tmp[[7]][["anc"]]]
# tmpy[is.na(tmpy)]<-cont.roots
# plot(0,type="n",xlim=range(node.depth.edgelength(tmp.tree)),ylim=range(seed))
# segments(x0=node.depth.edgelength(tmp.tree)[tmp.tree$edge[,1]],
#          x1=node.depth.edgelength(tmp.tree)[tmp.tree$edge[,2]],
#          y0=tmpy,y1=seed,
#          col=c("red","blue","green")[apply(maps,1,which.max)])
#looks like it works!

#so one potentially annoying quirk is that it needs multiple trees to do multiple discrete character
#histories...but that's as simple as adding some kind of on/off switch to assume all different or all
#same histories I think
# tree<-phytools::pbtree(n=100,scale=1)
# test<-sim.sce(rep(tree,100),nstates=2,freqs=c(0,5),jumps=9,nsims=100)
# phytools::phenogram(attr(test,"param_info")$params$trees[[2]],
#                     test$cont[,2],
#                     ftype="off")

#function to get parameter info from sce_sim objects
#works fine for now
#' @export
get.params<-function(sce.sim,
                     params=c("trees",
                              "disc.mods","cont.mods",
                              "rates","drifts",
                              "freqs","jumps","skews",
                              "disc.roots","cont.roots"),
                     sim=NULL){
  if(!inherits(sce.sim,"sce_sim")){
    stop("sce.sim must be an object of class 'sce_sim'")
  }

  params.choices<-c("trees",
                    "disc.mods","cont.mods",
                    "rates","drifts",
                    "freqs","jumps","skews",
                    "disc.roots","cont.roots")
  if(is.character(params)){
    params<-pmatch(params,params.choices)
  }else if(is.numeric(dist)){
    params[params>length(params.choices)]<-NA
  }else{
    stop("Didn't recognize parameter choices")
  }
  #should probably make a warning about this...
  params<-params[!is.na(params)]

  info<-attr(sce.sim,"param_info")
  nsims<-nrow(info[["lookups"]])
  if(is.null(sim)){
    sim<-seq_len(nsims)
  }else if(!is.numeric(sim)){
    stop("Didn't recognize sim; should be a vector of numeric indices")
  }
  #need a warning for this...
  sim<-sim[sim<=nsims]

  setNames(
    lapply(params.choices[params],
           function(ii) switch(ii,
                               disc.mods=info[["params"]][[ii]][,,info[["lookups"]][sim,ii],drop=FALSE],
                               cont.roots=info[["params"]][[ii]][info[["lookups"]][sim,ii]],
                               trees=info[["params"]][[ii]][info[["lookups"]][sim,ii]],
                               info[["params"]][[ii]][,info[["lookups"]][sim,ii],drop=FALSE])),
    params.choices[params])

  # if(length(out)==1){
  #   out[[1]]
  # }else{
  #   out
  # }
}
