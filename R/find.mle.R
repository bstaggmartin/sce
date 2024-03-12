#temporary maximum likelihood search function
#' @export
.fix.bds<-function(bds,upper,par.nms){
  if(is.list(bds)) bds<-unlist(bds)
  out<-bds[par.nms]
  nms<-names(bds)
  probs<-is.na(out)
  if(any(probs)){
    if(is.null(nms)){
      n.bds<-length(bds)
      nms<-character(n.bds)
      prob.nms<-!logical(n.bds)
    }else{
      nms[is.na(nms)]<-""
      prob.nms<-!(nms%in%par.nms)
    }
    if(any(prob.nms)){
      out[probs]<-rep(bds[prob.nms],length.out=sum(probs))
    }else if(upper){
      out[probs]<-Inf
    }else{
      out[probs]<- -Inf
    }
  }
  if(upper){
    out[is.infinite(out)&out<0]<-Inf
  }else{
    out[is.infinite(out)&out>0]<- -Inf
  }
  out
}

.get.init.fxn<-function(lb,ub,npar,init.width){
  #now forces all parameters to initialize around 0 as much as possible
  if(is.null(lb)) lb<-rep(-Inf,npar)
  if(is.null(ub)) ub<-rep(Inf,npar)
  lb.inf<-is.infinite(lb)
  ub.inf<-is.infinite(ub)
  lb[lb.inf]<-pmin(-init.width/2,ub[lb.inf]-init.width)
  ub[ub.inf]<-pmax(init.width/2,lb[ub.inf]+init.width)
  cents<-cbind((lb+ub)/2,lb+init.width/2,ub-init.width/2)
  whichs<-apply(abs(cents),1,which.min)
  cents<-cents[cbind(seq_len(npar),whichs)]
  lb<-pmax(lb,cents-init.width/2)
  ub<-pmin(ub,cents+init.width/2)
  function(inds){
    runif(length(inds),lb[inds],ub[inds])
  }
}

#allow for sequences of NLOPT calls by passing vectors to opts...
#just copied and minimally modified from contsimmap
#might want to update some things to make it a bit more intuitive going forward...
#' @export
find.mle<-function(lik.fun,init=NULL,times=1,lb=NULL,ub=NULL,...,
                   recycle.init=FALSE,init.width=5,
                   on.fail=c("random.restart","NA"),max.tries=100,
                   verbose=FALSE,step=1e-4){
  par.nms<-as.character(seq_len(max(lik.fun$param_key)))
  npar<-length(par.nms)
  if(!is.null(lb)){
    lb<-.fix.bds(lb,upper=FALSE,par.nms)
  }
  if(!is.null(ub)){
    ub<-.fix.bds(ub,upper=TRUE,par.nms)
  }
  #swap any bounds needing it
  if(!is.null(lb)&!is.null(ub)){
    tmp.lb<-lb
    tmp.ub<-ub
    tmp<-tmp.lb>tmp.ub
    lb[tmp]<-tmp.ub[tmp]
    tmp<-tmp.ub<tmp.lb
    ub[tmp]<-tmp.lb[tmp]
  }
  out.init<-matrix(nrow=npar,ncol=times,
                   dimnames=list(par.nms,NULL))
  if(!is.null(init)){
    if(!is.list(init)){
      ndims<-length(dim(init))
      if(ndims<2){
        if(is.null(names(init))){
          init<-list(init)
        }else{
          init<-split(init,names(init))
        }
      }else if(ndims>2){
        stop("Please input init as either a vector, matrix (with columns corresponding to different parameters), or list")
      }else{
        init<-asplit(init,2)
      }
    }
    nms<-names(init)
    if(is.null(nms)){
      n.init<-length(init)
      nms<-character(n.init)
      prob.nms<-!logical(n.init)
    }else{
      nms[is.na(nms)]<-""
      prob.nms<-!(nms%in%par.nms)
    }
    for(i in nms[!prob.nms]){
      if(recycle.init){
        out.init[i,]<-rep(init[[i]],length.out=times)
      }else{
        tmp.seq<-seq_len(min(times,length(init[[i]])))
        out.init[i,tmp.seq]<-init[[i]][tmp.seq]
      }
    }
    probs<-which(!(par.nms%in%nms[!prob.nms]))
    if(any(probs)){
      counter<-0
      if(recycle.init){
        for(i in rep(which(prob.nms),length.out=length(probs))){
          counter<-counter+1
          out.init[probs[counter],]<-rep(init[[i]],length.out=times)
        }
      }else{
        for(i in which(prob.nms)[seq_len(min(sum(prob.nms),length(probs)))]){
          counter<-counter+1
          tmp.seq<-seq_len(min(times,length(init[[i]])))
          out.init[probs[counter],tmp.seq]<-init[[i]][tmp.seq]
        }
      }
    }
  }
  init<-out.init
  nas<-is.na(init)
  inds<-(which(nas)-1)%%npar+1
  init.fxn<-.get.init.fxn(lb,ub,npar,init.width)
  init[nas]<-init.fxn(inds)
  opts<-list(...)
  if(is.null(opts[["algorithm"]])){
    opts[["algorithm"]]<-"NLOPT_LD_LBFGS"
  }
  if(is.null(opts[["maxeval"]])){
    opts[["maxeval"]]<-1e4
  }
  if(is.null(opts[["ftol_rel"]])){
    opts[["ftol_rel"]]<-sqrt(.Machine$double.eps)
  }
  if(is.null(opts[["xtol_res"]])){
    opts[["ftol_rel"]]<-sqrt(.Machine$double.eps)
  }
  if(length(on.fail)==1){
    if(is.na(on.fail)) on.fail<-"NA"
  }
  if(is.character(on.fail)) on.fail<-pmatch(on.fail[1],c("random.restart","NA"))
  if(is.na(on.fail)|!is.numeric(on.fail)) on.fail<-1
  if(verbose&is.null(opts[["print_level"]])){
    opts[["print_level"]]<-3
  }
  runs.per.time<-max(lengths(opts))
  opts<-lapply(opts,function(ii) rep(ii,length.out=runs.per.time))
  pars<-matrix(nrow=npar,ncol=times)
  codes<-liks<-numeric(times)
  for(i in seq_len(times)){
    if(verbose) cat("Maximizing likelihood function... (rep ",i," out of ",times,"):\n\n",sep="")
    for(j in seq_len(runs.per.time)){
      if(verbose) cat("Optimization round ",j," out of ",runs.per.time,"... (rep ",i," out of ",times,"):\n\n",sep="")
      tmp.opts<-lapply(opts,'[[',j)
      if(grepl("_LD_",tmp.opts[["algorithm"]])){
        res<-nloptr::nloptr(init[,i],lik.fun[["grad_lik"]],lb=lb,ub=ub,opts=tmp.opts,step=step)
      }else{
        res<-nloptr::nloptr(init[,i],lik.fun[["lik"]],lb=lb,ub=ub,opts=tmp.opts)
      }
      init[,i]<-res[["solution"]]
    }
    if(on.fail==1){
      counter<-1
      while(counter<=max.tries&(is.infinite(res[["objective"]])|res[["status"]]<0)){
        if(verbose) cat("Random restart ",counter," out of ",max.tries,"... (rep ",i," out of ",times,"):\n\n",sep="")
        counter<-counter+1
        init[,i]<-init.fxn(seq_len(npar))
        for(j in seq_len(runs.per.time)){
          if(verbose) cat("Optimization round ",j," out of ",runs.per.time,"... (rep ",i," out of ",times,"):\n\n",sep="")
          tmp.opts<-lapply(opts,'[[',j)
          if(grepl("_LD_",tmp.opts[["algorithm"]])){
            res<-nloptr::nloptr(init[,i],lik.fun[["grad_lik"]],lb=lb,ub=ub,opts=tmp.opts,step=step)
          }else{
            res<-nloptr::nloptr(init[,i],lik.fun[["lik"]],lb=lb,ub=ub,opts=tmp.opts)
          }
          init[,i]<-res[["solution"]]
        }
      }
    }
    liks[i]<- -res[["objective"]]
    pars[,i]<-res[["solution"]]
    codes[i]<-res[["status"]]
  }
  tmp.inds<-which(!(codes<0))
  #should also check for max iteration warnings
  #should rounding errors (-4) also be reported separately?
  if(!length(tmp.inds)){
    warning("Optimization algorithm failed to properly converge")
    final.max<-which.max(liks)
  }else{
    final.max<-tmp.inds[which.max(liks[tmp.inds])]
  }

  list("lnLik"=liks[final.max],
       "estimates"=pars[,final.max],
       "return_code"=codes[final.max])
}
