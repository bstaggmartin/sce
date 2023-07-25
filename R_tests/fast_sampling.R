#assumes size is 1 and n is equal to dim(prob)[1]
vec.sample<-function(prob){
  d<-NROW(prob)
  n<-NCOL(prob)
  cumprob<-apply(matrix(prob,d,n),2,cumsum)
  seed<-runif(n,max=cumprob[d,])
  foo<-function(ii){
    findInterval(seed[ii],cumprob[,ii])
  }
  out<-unlist(lapply(seq_len(n),foo),use.names=FALSE)
  out
}
tmp<-matrix(rexp(100*1024),1024,100)
vec.sample(tmp)
microbenchmark::microbenchmark(vec.sample(tmp))
microbenchmark::microbenchmark(unlist(lapply(seq_len(100),function(ii) sample(1024,1,prob=tmp[ii,])),use.names=FALSE))
#vec.sample cuts time by ~half
vec.sample2<-function(prob){
  d<-NROW(prob)
  n<-NCOL(prob)
  sums<-.colSums(prob,d,n)
  prob<-sweep(prob,2,sums,'/')
  out<-vector(length=n)
  i<-0
  seed<-runif(n)
  tmp.prob<-rep(0,n)
  avail.inds<-seq_len(n)
  while(length(avail.inds)){
    i<-i+1
    tmp.prob<-tmp.prob+prob[i,avail.inds]
    nopes<-seed>tmp.prob
    if(any(!nopes)){
      out[avail.inds[!nopes]]<-i
      seed<-seed[nopes]
      tmp.prob<-tmp.prob[nopes]
      avail.inds<-avail.inds[nopes]
    }
  }
  out
}
vec.sample2(tmp)
microbenchmark::microbenchmark(vec.sample2(tmp)) #slower than vec.sample, faster than sample
vec.sample3<-function(prob){
  d<-NROW(prob)
  n<-NCOL(prob)
  sums<-.colSums(prob,d,n)
  cumprob<-apply(matrix(prob,d,n),2,cumsum)
  seed<-runif(n,max=cumprob[d,])
  out<-vector(length=n)
  i<-0
  avail.inds<-seq_len(n)
  while(length(avail.inds)){
    i<-i+1
    nopes<-seed>cumprob[i,avail.inds]
    if(any(!nopes)){
      out[avail.inds[!nopes]]<-i
      seed<-seed[nopes]
      avail.inds<-avail.inds[nopes]
    }
  }
  out
}
vec.sample3(tmp)
microbenchmark::microbenchmark(vec.sample3(tmp)) #about the same speed as vec.sample2

tmp[sample(length(tmp),10000)]<-0
vec.sample(tmp) #seems to work with categories of 0...

#further testing indicated some numerical underflow issues with my "fast" approach


prob<-matrix(rexp(1024*100),1024,100)
prob[sample(length(prob),100000)]<-0
fast.sample(prob)
microbenchmark::microbenchmark(fast.sample(prob)) #about the same as sample...but I guess cleaner
microbenchmark::microbenchmark(unlist(lapply(seq_len(100),function(ii) sample(1024,1,prob=prob[,ii])),use.names=FALSE))
base.sample<-function(prob){
  unlist(lapply(seq_len(100),function(ii) sample(1024,1,prob=prob[,ii])),use.names=FALSE)
}
res.sample<-sapply(seq_len(1000),function(ii) base.sample(prob))
hist(res.sample[2,],breaks=200)
res.sample2<-sapply(seq_len(1000),function(ii) fast.sample(prob))
hist(res.sample2[2,],breaks=200)
plot(tabulate(res.sample2[2,],nbins=1024)~prob[,2])
plot(tabulate(res.sample[2,],nbins=1024)~prob[,2]) #both seem to work
#Based on how time varies with sparsity when there's a lot of 0s, I think sample() is doing something similar to what I'm doing "under the hood"
