library(rexpokit)
library(expm)
Q<-rand.Q(k,5,5)
diag(Q)<-diag(Q)-c(10,5)
t<-100
microbenchmark::microbenchmark(expm::expm(t*Q)%*%c(0.123,1-0.123))
microbenchmark::microbenchmark(rexpokit::expm(t*Q)%*%c(0.123,1-0.123))
eig<-eigen(Q)
eigval<-eig$values
eigvec<-eig$vectors
inveigvec<-solve(eigvec)
microbenchmark::microbenchmark((sweep(eigvec,2,exp(t*eigval),'*',check.margin=FALSE)%*%inveigvec)%*%c(0.123,1-0.123)) #fastest
microbenchmark::microbenchmark((sweep(eigvec,2,exp(t*eigval),'*',check.margin=FALSE))%*%(inveigvec%*%c(0.123,1-0.123)))


#too inaccurate, though the expm function is a bit faster using rexpokit!
# microbenchmark::microbenchmark(rexpokit::dgexpv(t(Q),input_probs_c(0.123,1-0.123))
# cooQ<-mat2coo(t(Q))
# debug(rexpokit::expokit_dgexpv_Qmat)
# microbenchmark::microbenchmark(rexpokit::expokit_dgexpv_Qmat(cooQ,t=1,inputprobs_for_fast=c(0.123,1-0.123),
#                                                              transform_to_coo_TF=FALSE,coo_n=2))
#I honestly think the diagonalization approach is more generalizable and just as numerically stable...

tmp<-XX[des[[i]],,,drop=FALSE]
nx<-2048
x<-apply(tmp,c(1,2),function(ii) IFFT(ii),simplify=FALSE)
holder2D<-matrix(0,2,nx)
# x<-apply(tmp,c(1,2),IDFT,simplify=FALSE)
for(i in seq_len(ncol(x))){
  holder2D[i,]<-Reduce('*',x[,i])
}
plot(abs(holder2D[1,]),type='l')
plot(abs(holder2D[2,]),type='l')
test<-t(apply(holder2D,1,FFT))
plot(abs(test[1,]),type='l')
plot(abs(test[2,]),type='l')
IDFT.test<-t(apply(test,1,IFFT))

#"manual" convolution --> only seems slightly more numerically stable, but officially several orders of magnitude faster :P
test.out<-tmp[1,,]
test.out[]<-0
x.seq<-seq_len(nx)
seqs<-lapply(x.seq,function(ii) ((ii-x.seq)%%nx)+1)
#non-circular convolution = bad idea, apparently
#so then why the difference? And is one better than the other?
# foo<-function(x){
#   out<-x-x.seq
#   out[out<=0|out>nx]<-NA
#   out
# }
# seqs<-lapply(x.seq,foo)
# ff<-tmp[1,,]
# gg<-tmp[2,,]
for(i in x.seq){
  test.out[,i]<-rowMeans(ff*gg[,seqs[[i]],drop=FALSE])
}
plot(abs(test.out[1,]),type='l')
plot(abs(test.out[2,]),type='l')
IDFT.test.out<-apply(test.out,1,IFFT)
plot(abs(IDFT.test.out[,1]),type='l')
plot(abs(IDFT.test.out[,2]),type='l')

matplot(log(abs(cbind(t(IDFT.test)))),type='l')
matplot(log(abs(IDFT.test.out)),type='l')

matplot(Re(t(test.out)),type='l')
matplot(Im(t(test.out)),type='l')

matplot(log(abs(t(holder2D))),type='l')
matplot(log(abs(IDFT.test.out)),type='l')

#there is a vectorized version, but looking at this example (which is still a lot slower than FFT convolution)...
#I'm not sure the improved numerical stability is universally true
#Though I would like to figure out where the difference in scaling comes from!
#maybe you're not technically doing circular convolution? (which would probs be better anyways?)
holder<-array(0,c(2,nx,nx))
tmptmp<-unlist(seqs,use.names=FALSE)
holder[1,,]<-ff[1,]*gg[1,tmptmp]
holder[2,,]<-ff[2,]*gg[2,tmptmp]
test.out2<-apply(holder,c(1,3),sum)
plot(abs(test.out2[1,]),type='l')
plot(abs(test.out2[2,]),type='l')

plot(as.vector(abs(IDFT.test.out/t(IDFT.test))),ylim=c(0,10000))
abline(h=2048)
#interesting...seems like manual convolution just scaled up by 2048. Factor must compound on itself somehow...
#fixed by taking rowMeans, rather than rowSums
#I would say that any potential numerical advantage of direct convolution in frequency domain is not worth the additional computation
#though, comparing the frequency domain functions specifically, there does seem to definitely be some advantage:
matplot(log(cbind(t(abs(test)),t(abs(test.out)))),type='l') #notable "truncation" of convolved distributions
#essentially means that a lot of error comes from FFT itself!
#but would slow things down by a factor of 400 for nx=1024, not to mention lack of zero padding/cutting things off at each convolution step...
#but at least now you know!

#one last try--matrix multiplication?
holder2<-matrix(gg[1,tmptmp],nx,nx)
holder3<-matrix(0,2,nx)
holder3[1,]<-ff[1,]%*%holder2
holder2<-matrix(gg[2,tmptmp],nx,nx)
holder3[2,]<-ff[2,]%*%holder2
holder3<-holder3/2048
plot(abs(holder3[1,]),type='l')
matplot(log(t(abs(test.out2))),type='l')
matplot(log(t(abs(holder3))),type='l')
matplot(log(t(abs(test.out))),type='l')
microbenchmark::microbenchmark({holder2<-matrix(gg[1,tmptmp],nx,nx)
holder3<-matrix(0,2,nx)
holder3[1,]<-ff[1,]%*%holder2
holder2<-matrix(gg[2,tmptmp],nx,nx)
holder3[2,]<-ff[2,]%*%holder2
holder3<-holder3/2048})
#cuts down to about 200x slower than FFT...better, but not quite there
holder3<-matrix(0,2,nx)
microbenchmark::microbenchmark({holder2<-array((gg/2048)[,tmptmp,drop=FALSE],c(2,nx,nx))
holder3[1,]<-ff[1,,drop=FALSE]%*%holder2[1,,]
holder3[2,]<-ff[2,,drop=FALSE]%*%holder2[2,,]},
times=10) #this is actually slower somehow...in between vectorized and above matrix mult strategy--maybe something to do with array construction and indexing?
#and you would have to repeat for each descendant edge...ugh

#SUMMARY: direct convolution in frequency domain is more numerically stable (due to underflow errors regarding FFT--seems to truncate any values
##lower than 1e-17th of the max value [about -40 on natural log scale])
#On the other hand, direct convolution is, at best, 200 times slower than FFT-based convolution
#Right now, I would say the cost of direct convolution via matrix multiplication outweighs the benefits since this only affects the tails of
##the distributions (exp(-40) or 1e-17 is a VERY small number!)
#It's also worth mentioning that this truncation threshold might vary according to a bunch of different factors...
#Warrants further exploration, I suppose
#On the other hand, this does create a slightly insidious side effect where discrete states with higher peaks will have higher tail densities
##since the truncation threshold varies with the peak
#You might get at this by actually re-scaling densities in each discrete state separately such that they both peak at 1...
#But this effect is probably very minor in practice!
#Still, probably worth exploring keeping track of separate scalars for each discrete state...than just calculate marginal likelihoods for each
##discrete separately, THAN marginalize over the resulting k likelihoods
#Probably would have relatively little speed impact but improve numerical stability in the tails!
