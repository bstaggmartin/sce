#3/4 notes:
# - gradient-based optimization is decently fast, but doesn't seem to offer speed benefits over subplex-based
# - putting bounds on parameter estimates seems to lead to headaches more often than not
# - I think a default init.width of 5 is fine given the log scale of most parameters...
# --> definitely seems to imporve fitting behavior, though I worry this could break down in certain cases...

#3/5 notes:
# - using a simple forward finite diff approx is fast and pretty accurate here
# - seems better tha subplex...
# - init.width of 5 seems like a pretty good idea (note that exp(10) spans range of 22,000)
# - no bounds necessary I think
# - the only issue is unidentifiability in some parts of param space when Levy processes enter the picture
# - if freq is close to 0 jump has virtually no effect on likelihood and vice versa
# - how to fix that?
# - could levy processes be characterized by their overall variance and an additional parameter?
# - probs not that simple, but worth thinking on some more...
# - yeah, definitely some weak identifiability issues at play here with mixing Levy and BM processes...but something to figure out in the future
# - skews might actually be weakly identifiable under certain conditions!
# - bobyqa algorithm also works pretty well!
# - BFGS definitely the most consistent

#3/7 notes:
# - NIG simulation now fully working again
# - interestingly, the freq vs jump parameter of NIG processes seem somewhat unidentifiable!
# - like, I think their product is identifiable, but the specific numbers are difficult to pin down

#3/11 notes:
# - both vg and nig process simulation and fitting all synced up!
# - vg and nig skews very weakly identifiable
# - I made the delta param equal to freq*jump for greater consistency/interp
# ...but they're still only weakly identifiable
# - probably make an alpha stable process next

library(sce)

set.seed(321)
tree<-phytools::pbtree(n=100,scale=1)
test<-sim.sce(tree,
              nstates=4,
              disc.mods=matrix(c(-2,1,1,0,1,-2,0,1,1,0,-2,1,0,1,1,-2),4,4),
              cont.mods="JN",
              rates=c(1,5,5,10),
              drifts=c(0,0),
              freqs=c(0,0),
              jumps=c(0,0),
              skews=c(0,0),
              disc.roots=1)
phytools::phenogram(get.params(test,"tree")[[1]][[1]],test$cont[,1],ftype="off")
cont<-test$cont[,1]
disc<-test$disc[,1]

new.disc<-setNames(
  rep(c("STATE_1A&STATE_1B","STATE_2A&STATE_2B"),each=2)[match(disc,c("STATE_1","STATE_2","STATE_3","STATE_4"))],
  names(disc)
)
# new.disc<-setNames(
#   do.call(paste,c(lapply(LETTERS[1:2],function(ii) paste0(disc,ii)),sep="&")),
#   names(disc)
# )


Qo<-matrix(c(0,1,2,0),2,2,
           dimnames=rep(list(sort(unique(disc))[1:2]),2))
Qh<-matrix(c(0,3,4,0),2,2,
           dimnames=rep(list(LETTERS[1:2]),2))
tmplt<-as.matrix(Matrix::bdiag(rep(list(Qh),nrow(Qo)))+kronecker(Qo,diag(nrow(Qh))))
rownames(tmplt)<-colnames(tmplt)<-paste0(rownames(Qo),rep(rownames(Qh),each=nrow(Qh)))
fun<-make.sce(tree,new.disc,cont,res=256,disc.mods=tmplt)

# tmp.seed<-rnorm(8)
# numDeriv::grad(fun$lik,tmp.seed)-fun$grad_lik(tmp.seed)$gradient
fit<-find.mle(fun,verbose=TRUE,lb=-10,ub=10)



fit<-find.mle(fun,verbose=TRUE,init.width=5,algorithm="NLOPT_LD_LBFGS")

find.mle(fun,init=t(fit$estimates),verbose=TRUE,init.width=5,algorithm="NLOPT_LD_LBFGS")

hess<-numDeriv::hessian(fun$lik,fit$estimates)
sqrt(diag(solve(hess)))
tmp.foo<-function(x) fun$grad_lik(x)[["gradient"]]
hess2<-numDeriv::jacobian(tmp.foo,fit$estimates,method="simple") #same but a lot quicker

tmp<-make.sce(tree,disc,cont,res=512,freqs="SD",jumps="SD")
test.fit<-nloptr::nloptr(rnorm(8),tmp[[1]],opts=list(algorithm="NLOPT_LN_SBPLX",
                                                     xtol_rel=sqrt(.Machine$double.eps),
                                                     ftol_rel=sqrt(.Machine$double.eps),
                                                     maxeval=1e6,
                                                     print_level=3))
# test.fit<-optim(rnorm(8),tmp[[1]],method="BFGS",control=list(trace=TRUE),hessian=TRUE)
# solve(test.fit$hessian)
rec<-fun$recon(fit$estimates,probs=c(0.5,0.025,0.975,0.001,0.999))

Cairo::CairoPDF("test_recon.pdf",width=8.5,height=11)
phytools::phenogram(#tree,
  get.params(test,"tree")[[1]][[1]],
  colors=setNames(hcl.colors(5)[-5],paste0("STATE_",c(2,1,4,3))),
  rec$cont[,1],ftype="off",
  ylim=range(rec$cont),
  lwd=2)
segments(x0=ape::node.depth.edgelength(tree),
         y0=rec$cont[,2],
         y1=rec$cont[,3])
segments(x0=ape::node.depth.edgelength(tree)-0.01,
         x1=ape::node.depth.edgelength(tree)+0.01,
         y0=c(rec$cont[,2],rec$cont[,3]))
ape::nodelabels(pie=rec$disc[-seq_len(ape::Ntip(tree)),],
                piecol=hcl.colors(5)[-5],
                cex=0.5)
ape::tiplabels(pie=rec$disc[seq_len(ape::Ntip(tree)),],
               piecol=hcl.colors(5)[-5],
               cex=0.5)
dev.off()

Cairo::CairoPDF("test_recon3.pdf",width=8.5,height=11)
#violin plot test...
phytools::phenogram(#tree,
                    get.params(test,"tree")[[1]][[1]],
                    colors=setNames(hcl.colors(5)[-5],paste0("STATE_",c(2,1,4,3))),
                    rec$cont[,1],ftype="off",
                    ylim=range(rec$cont),
                    lwd=2)
wd<-diff(par()$usr[1:2])
xx<-attr(rec$joint,"xpts")
cols<-hcl.colors(5)[-5]
for(i in seq_len(dim(rec$joint)[1])){
  tmp.inds<-xx>rec$cont[i,4]&xx<rec$cont[i,5]
  # tmp.inds<-rep(TRUE,dim(rec$joint)[3])
  tmp.yy<-xx[tmp.inds]
  tmp<-rec$joint[i,,tmp.inds]
  marg.tmp<-colSums(tmp)
  tmp<-tmp/max(marg.tmp)*wd/50
  marg.tmp<-marg.tmp/max(marg.tmp)*wd/50
  xx0<-tmp.xx<-ape::node.depth.edgelength(tree)[i]-marg.tmp/2
  for(j in seq_len(dim(rec$joint)[2])){
    new.xx<-tmp.xx+tmp[j,]
    polygon(c(tmp.yy,rev(tmp.yy))~c(tmp.xx,rev(new.xx)),
            col=cols[j],border=NA)
    tmp.xx<-new.xx
  }
  polygon(c(tmp.yy,rev(tmp.yy))~c(xx0,rev(tmp.xx)),density=0)
}
dev.off()


matplot(t(rec$dists[111,,]),x=attr(rec$dists,"xpts"),type="l")
abline(v=rec$dists_sum[111,4:5])

tmpy<-phytools::contMap(tree,rec$dists_sum[seq_len(ape::Ntip(tree)),3],
                        anc.states=rec$dists_sum[-seq_len(ape::Ntip(tree)),3])
tmpy<-phytools::setMap(tmpy,hcl.colors(10))
plot(tmpy,ftype="off")
ape::nodelabels(pie=rec$dists_sum[-seq_len(ape::Ntip(tree)),1:2],
                piecol=hcl.colors(2,palette="Dark3"),
                cex=0.5)
ape::tiplabels(pie=rec$dists_sum[seq_len(ape::Ntip(tree)),1:2],
               piecol=hcl.colors(2,palette="Dark3"),
               cex=0.5)

plot(get.params(test,"trees")[[1]][[1]])
ape::nodelabels(pie=rec$dists_sum[-seq_len(ape::Ntip(tree)),1:2],
                piecol=1:2,
                cex=0.5)
ape::tiplabels(pie=rec$dists_sum[seq_len(ape::Ntip(tree)),1:2],
               piecol=1:2,
               cex=0.5)

#state probs
marg.probs<-apply(test,c(1,2),sum)
#means
# state.means<-apply(test,c(1,2),function(ii) sum(ii*tmp[[3]])/sum(ii))
# #"unavailable" means coded as NaNs!
# marg.means<-rowSums(marg.probs*state.means,na.rm=TRUE)

#simpler with no NaNs:
marg.means<-apply(test,1,function(ii) sum(tmp[[3]]*colSums(ii)))

#confidence intervals
ii<-test[1,,]
dx<-tmp[[3]][2]-tmp[[3]][1]
foo<-function(ii){
  cdf<-cumsum(colSums(ii))
  lb<-min(which(cdf>0.025))
  ub<-min(which(cdf>0.975))
  c((0.975-cdf[ub-1])*dx/(cdf[ub]-cdf[ub-1])+tmp[[3]][ub-1],
    (0.025-cdf[lb-1])*dx/(cdf[lb]-cdf[lb-1])+tmp[[3]][lb-1])
}
apply(test,1,foo)

matplot(log(t(test[,2,])),type="l")

plot(cumsum(test[1,2,]))


#apparent variation in lik as resolution/bounds change is due to cont.se getting rounded up to 3*dx only!
#otherwise, log likelihoods are generally identical to 2-4 decimal places
#(at least for 2-rate BM...)
#yeah, decided to stop adjusting cont.se based on resolution and instead just return warnings!
#it's not great for Levy models it seems (liks change pretty drastically depending on res based on
#effective tip error)...but better/more intuitive than auto-adjusting
#standard errors now default to 1/100th of trait interval
#seems to work pretty well
#even with drifts and such!
likfun<-make.sce(tree,disc,cont,res=256)
likfun(log(c(1,1,1,9)))
test.fit<-optim(rnorm(4),likfun,method="BFGS",hessian=TRUE)
exp(sqrt(diag((solve(test.fit$hessian)))))

#from an old test with drifts:
#wideeee intervals on drift parameter as you'd expect!
#and yeah, seems like relative differences in drift are weakly identifiable, but not overall drift!

library(dentist)
undebug(dent_walk)
test.ci<-dent_walk(test.fit$par,likfun,test.fit$value,lower_bound=-Inf,nsteps=5000)
exp(test.ci$all_ranges)
plot(test.ci$results)
plot(test.ci$results[,c(4,5)])
