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

sim.trajectory<-function(t,x0,xlim,sig2,theta,alpha,lambda,gam2,res=100){
  # xlim<-xlim[!is.infinite(xlim)&!is.na(xlim)]
  # if(length(xlim)<2){
  #   stop("xlim must consist of 2 numbers")
  # }
  # xlim<-sort(xlim[seq_len(2)])
  # if(x0<xlim[1]|x0>xlim[2]){
  #   stop("x0 must fall with xlim interval")
  # }
  tpts<-seq(0,t,length.out=res+1)
  jumps<-NULL
  tt<-rexp(1,lambda)
  while(tt<t){
    jumps<-c(jumps,tt)
    tt<-tt+rexp(1,lambda)
  }
  tpts<-c(jumps,tpts)
  jumps<-rep(c(TRUE,FALSE),c(length(jumps),res+1))
  ord<-order(tpts)
  tpts<-tpts[ord]
  jumps<-jumps[ord]
  dts<-diff(tpts)
  nt<-length(tpts)
  out<-vector('numeric',length=nt)
  out[1]<-x0
  seed<-rnorm(nt-1+sum(jumps))
  counter<-1
  for(i in seq_along(dts)){
    tmp<-exp(-alpha*dts[i])
    out[i+1]<-sqrt(sig2/(2*alpha)*(1-tmp^2))*seed[counter]+out[i]*tmp+theta*(1-tmp)
    counter<-counter+1
    if(jumps[i+1]){
      out[i+1]<-out[i+1]+sqrt(gam2)*seed[counter]
      counter<-counter+1
    }
    tmp<-c(out[i+1]<xlim[1],out[i+1]>xlim[2])
    while(any(tmp)){
      out[i+1]<-2*xlim[tmp]-out[i+1]
      tmp<-c(out[i+1]<xlim[1],out[i+1]>xlim[2])
    }
  }
  list(out[nt],sum(dts*exp(out[-1])))
}

sim.evorates<-function(tree,x0,r0,rlim,rsig2,rtheta,ralpha,rlambda,rgam2,res=100){
  list2env(sce:::.get.tree.topo.info(tree,pruning=FALSE),envir=environment())

  ##discretizing time##
  #get tree height and edge lengths
  hgt<-max(node.depth.edgelength(tree))
  #"ideal" time interval based on nt
  dt<-hgt/res
  #find time res along each branch
  nts<-ceiling(elen/dt)

  nnode<-max(nodes,na.rm=TRUE)
  ntip<-length(tips)
  rr<-xx<-vector('numeric',nnode)
  xx[ntip+1]<-x0
  rr[ntip+1]<-r0
  bw.rates<-vector('numeric',nedge-1)
  for(i in clade.seq){
    tmp.nodes<-nodes[i,]
    x0<-xx[tmp.nodes[1]]
    r0<-rr[tmp.nodes[1]]
    tmp<-sim.trajectory(elen[i],r0,rlim,rsig2,rtheta,ralpha,rlambda,rgam2,res=nts[i])
    rr[tmp.nodes[2]]<-tmp[[1]]
    bw.rates[i]<-tmp[[2]]
    xx[tmp.nodes[2]]<-x0+rnorm(1,0,sqrt(bw.rates[i]))
  }

  nms<-c(tree[['tip.label']],seq_len(nnode-ntip)+ntip)
  list('tree'=tree,
       'branchwise.rates'=bw.rates,
       'X'=setNames(xx,nms),
       'R'=setNames(rr,nms))
}

library(phytools)
tree<-pbtree(n=75,scale=1)
test<-sim.evorates(tree,x0=0,r0=0,rlim=c(-10,10),rsig2=1,rtheta=-4,ralpha=1,rlambda=10,rgam2=2)
phenogram(test[['tree']],test[['X']],ftype='off')

#seems to workish...still some kinks to iron out by the looks of it
#in particular, gam2 seems to influence the jump frequency for some reason...maybe because I standardize scaling?
#nope, just a problem with the normal DFT--now seems to work!
#I think your definition of directionality works, but it seems to lead to aliasing issues! Drat.
#Seems to necessitate a cleaning step to get rid of "off-peak" fluctuations...
#wayyyy better when you directly exponentiate the matrix! Also mu seems to have to be divided by 2 for some reason?
#having a mu for jumps is still an unsolved issue--would have to compute the jump kernels in a drastically different way
nx<-512
x0<- -10
x1<-20
dx<-(x1-x0)/(nx-1)
alpha<-0.4
theta<-5
sig2<-1
lambda<-1
gam2<-100
get.Q<-get.Q.fxn(nx)
Q<-get.Q(x0,dx,theta,alpha,sig2,lambda,gam2)
# get.P<-get.Q(dx,mu,sig2,lambda,gam2)
dt<-0.05
t<-50
nt<-t/dt+1
P<-expm::expm(dt*Q)
P[P<0]<-0

x<-vector('integer',nt)
x[1]<-sample(nx,1)
for(i in seq_len(nt)[-1]){
  x[i]<-sample(nx,1,prob=P[x[i-1],])
}
xx<-(x-1)*dx+x0
tt<-(seq_len(nt)-1)*dt
plot(xx~tt,type='l',ylim=c(x0,x1));abline(h=theta)

# #old
# get.Q.base<-function(nx){
#   base<-matrix(0,nx,nx)
#   rows<-.row(c(nx,nx))
#   cols<-.col(c(nx,nx))
#   diag.inds<-rows==cols
#   off.diag.inds<-rows==cols-1|cols==rows-1
#   xpts<-seq_len(nx)
#   shift.inds<-base
#   #doesn't quite work since you ideally want reflexive bounds, but could use FFT for this...
#   for(i in seq_len(nx)){
#     shift.inds[-seq_len(i),i]<-seq_len(nx-i)
#   }
#   shift.inds[rows<cols]<-t(shift.inds)[rows<cols]
#   function(dx,sig2,lambda,gam2){
#     Q2<-Q1<-base
#     Q1[off.diag.inds]<-sig2/(2*dx^2)
#     Q1[diag.inds]<- -.rowSums(Q1,nx,nx)
#     xpts<-dx*xpts
#     Q2[!diag.inds]<-diff(pnorm(xpts,0,sqrt(gam2)))[shift.inds]
#     Q2<-sweep(Q2,1,lambda/.rowSums(Q2,nx,nx),'*',check.margin=FALSE)
#     Q2[diag.inds]<- -lambda
#     eig<-eigen(Q1+Q2)
#     eigval<-eig[['values']]
#     eigvec<-eig[['vectors']]
#     inveigvec<-solve(eigvec)
#     function(t){
#       sweep(eigvec,2,exp(t*eigval),'*',check.margin=FALSE)%*%inveigvec
#     }
#   }
# }
#
# #older
# #maybe add some parameter such that jumps decrease in probability with distance?
# get.Q.base<-function(nx,dx){
#   base<-matrix(0,nx,nx)
#   rows<-.row(c(nx,nx))
#   cols<-.col(c(nx,nx))
#   diag.inds<-rows==cols
#   off.diag.inds<-rows==cols-1|cols==rows-1
#   K<-1/(2*dx^2)
#   function(sig2,lambda){
#     base[!diag.inds]<-lambda
#     base[off.diag.inds]<-lambda+K*sig2
#     base[diag.inds]<- -.rowSums(base,nx,nx)
#     eig<-eigen(base)
#     eigval<-eig[['values']]
#     eigvec<-eig[['vectors']]
#     inveigvec<-t(eigvec) #can only do this if symmetric, I believe
#     function(t){
#       sweep(eigvec,2,exp(t*eigval),'*',check.margin=FALSE)%*%inveigvec
#     }
#   }
# }
#
# #oldest
# BM.Q<-function(nx,dx){
#   base<-matrix(0,nx,nx)
#   inds<-seq_len(nx)
#   diag.inds<-cbind(inds,inds)
#   inds<-inds[-1]
#   inds2<-inds-1
#   off.diag.inds<-matrix(c(inds,inds2,inds2,inds),2*nx-2,2)
#   tmp<-1/(2*dx^2)
#   base[off.diag.inds]<-tmp
#   base[diag.inds]<- -2*tmp
#   base[c(1,nx^2)]<- -tmp
#   eig<-eigen(base)
#   eigval<-eig[['values']]
#   eigvec<-eig[['vectors']]
#   inveigvec<-t(eigvec)
#   out<-function(t,sig2){
#     sweep(eigvec,2,exp(t*sig2*eigval),'*',check.margin=FALSE)%*%inveigvec
#   }
# }
