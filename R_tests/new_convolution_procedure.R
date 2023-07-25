library(sce)

nx<-1024
dx<-0.1
xx<-seq(-nx/2*dx,by=dx,length.out=nx)
BM<-BM.DFT(nx,dx)
dirac<-dirac.DFT(nx,dx,x0=xx[1])
Q<-matrix(rexp(4),2,2)
diag(Q)<-0
diag(Q)<- -colSums(Q)
t<-1
nt<-101
dt<-t/(nt-1)
P<-expm::expm(dt*Q)
C<-cbind(BM(dt,1),BM(dt,100))
out1<-matrix(c(1,0),nrow=2*nx,ncol=ncol(Q),byrow=TRUE)
out2<-matrix(c(0,1),nrow=2*nx,ncol=ncol(Q),byrow=TRUE)
for(i in seq_len(nt-1)){
  out1<-C*(out1%*%P)
  out2<-C*(out2%*%P)
}
tmp.out1<-Re(apply(out1,2,fftw::IFFT))
tmp.out2<-Re(apply(out2,2,fftw::IFFT))
matplot(tmp.out1,type='l')
matplot(tmp.out2,type='l')
#now, I *think* any starting point is just a superposition of these base cases...right?

gauss<-gauss.DFT(nx,dx,x0=xx[1])
out<-init<-cbind(gauss(sample(xx,1),0.1),gauss(sample(xx,1),1))
for(i in seq_len(nt-1)){
  out<-C*(out%*%P)
}
out<-Re(apply(out,2,fftw::IFFT))
matplot(out,type='l')

#only issue is how to structure init...I think there might be something here, but it's not intuitive...
#you might actually just copy it over???
#oh my god yeah, it's actually the exact same!!!
#so you can come up with a convolution of sorts for an arbitrary starting point...holy shit
#just requires a different convolution for starting in any discrete state...
test<-Re(apply(init[,1]*out1+init[,2]*out2,2,fftw::IFFT))
matplot(test,type='l')
#so the idea, then, is to first construct a "lookup table"
#essentially approximate the convolutions for any discrete state along a sequence of time points, recording any relevent time points...
#would also be good for simply studing the structure of the convolutions...
#could this logic be extended beyond BM to things like OU? Proably not since it still relies on convolution at its core
#but neat connection to realize that bring further structure and clarity to the problem nonetheless!
#the cool thing is that it implies a theoretical characteristic function for any starting discrete state...
#but how to find that characteristic function remains a mystery...
#could always use matrix exponentiation, but that doesn't feel very helpful
#this only cuts down (admittedly SUBSTANTIALLY) on computation for your double-discretized framework
#think about it: you only have to convolve enough time points to get the longest branch in the phylogeny...
#and, wait, there's a recursion here...
#if you can form this into a convolution at an arbitrary time length, could you not just...use the current result to essentially "double" the last result?
#like, if you approximate the convolutions to dt, couldn't you get the 2*dt convolutions?
#Interesting...implies that a non-uniform time-discretization might be helpful...
#Or just a single, short-interval one that gets raised to a power and multiplied?
#But like a lot...is that just matrix exponentiation? lol

#basically would look like having...
out1<-C*matrix(P[,1],nrow=2*nx,ncol=ncol(Q),byrow=TRUE)
out2<-C*matrix(P[,2],nrow=2*nx,ncol=ncol(Q),byrow=TRUE)
#and somehow adding and multiplying these kernels to get any arbitrary time convolution...
#it feels like you may be able to take advantage of some kind of recursion to get that?
#so for getting further time points for state 1...
(out1[,1]*out1+out1[,2]*out2)[,1]*out[,1]+(out1[,1]*out1+out1[,2]*out2)[,2]*out2
#feels like dead end, but might be worth expanding out and seeing if there's any pattern to it
#I think the indexing will make it ineffectual...

#can this be simplified at all?
#not obviously...
new.out1<-out1[,1]*out1+out1[,2]*out2
# new.out[,1]<-out1[,1]*out1[,1]+out1[,2]*out2[,1]
# new.out[,2]<-out1[,1]*out1[,2]+out1[,2]*out2[,2]
#here's an interesting way of structuring it...
# new.out1[,1]<-rowSums(out1*cbind(out1[,1],out2[,1]))
# new.out1[,2]<-rowSums(out1*cbind(out1[,2],out2[,2]))
#implies keeping track of out1 and 2 as 3d array to be able to take slices across input OR output states

CC<-array(t(P),dim=c(ncol(Q),ncol(Q),2*nx))
CC<-sweep(CC,c(3,2),C,'*') #now dim1 corresponds to input state, dim2 to output state, and dim3 to continuous state
new.CC<-CC
new.CC[1,1,]<-colSums(CC[1,,]*CC[,1,])
new.CC[1,2,]<-colSums(CC[1,,]*CC[,2,])
new.CC[2,1,]<-colSums(CC[2,,]*CC[,1,])
new.CC[2,2,]<-colSums(CC[2,,]*CC[,2,])
#now some regularities are starting to become more apparent...
#I feel like you could turn this into a bunch of matrix multiplications for easier notation...
#hmmm...might be doing something wrong (probabilites don't seem to be conserved)
#but that might be a transposition thing?
#the really cool thing is that it looks like you're just taking matrix powers for each continuous state right now
new.CC[,,1];CC[,,1]%*%CC[,,1]
#only problem is that taking matrix powers isn't necessarilly easy either...
#or is it? If it's a small state space, would it be worth simply diagonalizing the caculations and using that?
#yeah, just a transposition thing-->I think I fixed it
#Interesting...so, if you start with this tiny little kernel, then diagonalize them for each state space, you end up with convolution
##matrices that can be combined to convolve for an arbitrary starting point and time interval...right?
tmp<-eigen(CC[,,1])
eigvec<-tmp$vectors
inveigvec<-solve(eigvec)
eigval<-tmp$values
#the below works
CC[,,1]%*%CC[,,1]
eigvec%*%diag(eigval^2)%*%inveigvec
#and this could DEFINITELY be vectorized for speedups--the main slow down would be diagnolizing all those matrices in the beginning!
#and it's unclear whether this could extend to non-integer powers? It looks like it can!
#This would be REALLY cool if it's not prone to numerical errors!



#welp, why not try it out? You've got time
holder<-CC<-sweep(array(P,dim=c(ncol(Q),ncol(Q),2*nx)),c(3,2),C,'*')
eigs<-apply(CC,3,eigen)
eigvec<-lapply(eigs,'[[','vectors')
inveigvec<-lapply(eigvec,solve)
holder[]<-unlist(eigvec,use.names=FALSE)
eigvec<-holder
holder[]<-unlist(inveigvec,use.names=FALSE)
inveigvec<-holder
holder[,1,]<-unlist(lapply(eigs,'[[','values'),use.names=FALSE)
eigval<-holder[,1,]
k<-ncol(Q)
#to get for arbitrary time point, t...
foo<-function(t){
  tmp<-sweep(eigvec,c(2,3),eigval^(t/dt),'*',check.margin=FALSE)
  #than figure out some way to matrix multiply each slice by inveigvec...
  #definitely some underflow stuff possible
  out<-tmp[,rep.int(1,k),,drop=FALSE]*inveigvec[rep.int(1,k),,,drop=FALSE]
  for(i in seq_len(k-1)){
    out<-out+tmp[,rep.int(i+1,k),,drop=FALSE]*inveigvec[rep.int(i+1,k),,,drop=FALSE]
  }
  out
}
test<-foo(1)
test[,,1]
Reduce('%*%',rep(list(CC[,,1]),10))
#I think I must just be doing the math wrong...
eigvec[,,1]%*%diag(eigval[,1]^10)%*%inveigvec[,,1]

# test[test<0]<- -test[test<0]
tmp.test1<-Re(apply(test[1,,],1,fftw::IFFT))
tmp.test2<-Re(apply(test[2,,],1,fftw::IFFT))
matplot(tmp.test1,type='l')
matplot(tmp.out1,type='l')
matplot(tmp.test2,type='l')
matplot(tmp.out2,type='l')

matplot(t(test[1,,]),type='l') #interesting-->some sort of numerical error from matrix diagonalization, I guess?
matplot(out1,type='l',add=TRUE)

matplot(t(test[2,,]),type='l') #interesting-->some sort of numerical error from matrix diagonalization, I guess?
matplot(out2,type='l',add=TRUE)
#only in second state it looks like
#oh...second eigen value gets "blown up" by the looks of it!
#tried flipping it and it didn't work! shoot...
#yeah, the diagonalization unfortunately doesn't seem numerically stable enough...
#maybe if you do it over a longer time period? Unsure...
#The issue probably has more to do with imperfect cancellation of eigenvectors since changing eigenvalues basically created no difference
#Unfortunately, increasing the timestep also yielded very little benefit
#There seems to be some issue where...
#Is it just a transposition issue?
#No that's not it...
#Just seems to be good old-fashioned numerical instability
#But I feel like there must be some way to rescue it just a bit?

#Is singular value decomp useful?
CC[,,1]
tmp<-svd(CC[,,1])
tmp$u%*%diag(tmp$d)%*%tmp$v
CC[,,1]%*%CC[,,1]
tmp$u%*%diag(tmp$d)^2%*%tmp$v

#so I figured out...something...it's close, but not quite there. Seems to be keeping things in state 1 for some reason...can't figure out why
#oh, nope! Just a transposition!
#this all works!!!
microbenchmark::microbenchmark(foo(100)) #takes a little over millisecond, regardless of time/resolution

#interesting observation:
#the characteristic function of a dirac delta centered at 0 is 1 everywhere
#so, a BM DFT essentially consists of a uniform distribution that gets expoentially "squashed" everywhere, with the middle getting
##squashed at a higher rate
#this clue gives some insight in how to derive a continuous time solution, but I think it really just reduces to the full matrix exponential
#Why? because the flux between discrete states is affected by the amount in each state, and the amount in each state is dynamically changing over
##in reponse to this squashing...
#But maybe you could encode it in the diagonals?
Q<-matrix(rexp(4),2,2)
diag(Q)<-0
diag(Q)<- -colSums(Q)
#the rate of squashing should be...
#base*sig2
#implicity something along the lines of sig2*x^2/2
base<-sce:::.get.base(1024,0.1)
Q
sig2<-c(1,100)
test<-array(Q,c(2,2,2048))-sweep(array(diag(sig2),c(2,2,2048)),3,base,'*')
testy<-test
testy[]<-apply(test,3,expm::expm)
tmp.test<-Re(apply(testy[1,,],1,fftw::IFFT))
matplot(tmp.test,type='l')
matplot(tmp.out1,type='l')
#hmmm...but how do I know which is "right"?
#both look similar, but they're definitely different!
big.Q<-matrix()


#formal test
#traditional way
nx<-1024
dx<-0.1
xx<-seq(-nx/2*dx,by=dx,length.out=nx)
BM<-BM.DFT(nx,dx)
dirac<-dirac.DFT(nx,dx,x0=xx[1])
Q<-matrix(rexp(4),2,2)
diag(Q)<-0
diag(Q)<- -colSums(Q)
t<-1
nt<-1001
dt<-t/(nt-1)
P<-expm::expm(dt*Q)
C<-cbind(BM(dt,1),BM(dt,100))
out1<-matrix(c(1,0),nrow=2*nx,ncol=ncol(Q),byrow=TRUE)
out2<-matrix(c(0,1),nrow=2*nx,ncol=ncol(Q),byrow=TRUE)
for(i in seq_len(nt-1)){
  out1<-C*(out1%*%P)
  out2<-C*(out2%*%P)
}
tmp.out1<-Re(apply(out1,2,fftw::IFFT))
tmp.out2<-Re(apply(out2,2,fftw::IFFT))
matplot(tmp.out1,type='l')
matplot(tmp.out2,type='l')
#new way
base<- -sce:::.get.base(nx,dx)/2
tmp<-seq(0,pi/dx,length.out=nx+1)
loc.base<- -1i*c(tmp,tmp[-c(1,nx+1)]+tmp[nx+1])
tmp<-tmp+diff(tmp[1:2])/2
loc.base2<- -1i*c(tmp,tmp[-c(1,nx+1)]+tmp[nx+1])
test<-array(Q,c(2,2,2*nx))+sweep(array(diag(c(0.1,3)),c(2,2,2*nx)),3,base,'*')+sweep(array(diag(c(1,-0.5)),c(2,2,2*nx)),3,loc.base,'*')
test2<-array(Q,c(2,2,2*nx))+sweep(array(diag(c(0.1,3)),c(2,2,2*nx)),3,base,'*')+sweep(array(diag(c(1,-0.5)),c(2,2,2*nx)),3,loc.base2,'*')
testy<-test
testy[]<-apply(test,3,expm::expm)
tmp.test<-Re(apply(testy[1,,],1,fftw::IFFT))
matplot(tmp.test,type='l')
matplot(tmp.out1,type='l')
#forgot a factor of 2--they're actually the same!!! Reduces whole problem to just nx k*k matrix exponentials...wow!
#super, super cool...
#so now you're only doing trait discretization, not time! I feel like folks will love this
#how fast, though?
microbenchmark::microbenchmark(apply(test,3,expm::expm))
#pretty slow --> matrix diagonalization could be numerically unstable, but beneficial here...
holder<-array(NA,c(2,2,2048))
eigs<-apply(test,3,eigen)
eigvec<-lapply(eigs,'[[','vectors')
inveigvec<-lapply(eigvec,solve)
holder[]<-unlist(eigvec,use.names=FALSE)
eigvec<-holder
holder[]<-unlist(inveigvec,use.names=FALSE)
inveigvec<-holder
holder[,1,]<-unlist(lapply(eigs,'[[','values'),use.names=FALSE)
eigval<-holder[,1,]
k<-ncol(Q)
#to get for arbitrary time point, t...
foo<-function(t){
  tmp<-sweep(eigvec,c(2,3),exp(t*eigval),'*',check.margin=FALSE)
  #than figure out some way to matrix multiply each slice by inveigvec...
  #definitely some underflow stuff possible
  out<-tmp[,rep.int(1,k),,drop=FALSE]*inveigvec[rep.int(1,k),,,drop=FALSE]
  for(i in seq_len(k-1)){
    out<-out+tmp[,rep.int(i+1,k),,drop=FALSE]*inveigvec[rep.int(i+1,k),,,drop=FALSE]
  }
  out
}
foo(1)

#looking into phase shift dealiasing...yeah, I just don't get the math...
eigs2<-apply(test2,3,eigen)
eigvec2<-lapply(eigs2,'[[','vectors')
inveigvec2<-lapply(eigvec2,solve)
holder[]<-unlist(eigvec2,use.names=FALSE)
eigvec2<-holder
holder[]<-unlist(inveigvec2,use.names=FALSE)
inveigvec2<-holder
holder[,1,]<-unlist(lapply(eigs2,'[[','values'),use.names=FALSE)
eigval2<-holder[,1,]
k<-ncol(Q)
#to get for arbitrary time point, t...
foo2<-function(t){
  tmp<-sweep(eigvec2,c(2,3),exp(t*eigval2),'*',check.margin=FALSE)
  #than figure out some way to matrix multiply each slice by inveigvec...
  #definitely some underflow stuff possible
  out<-tmp[,rep.int(1,k),,drop=FALSE]*inveigvec2[rep.int(1,k),,,drop=FALSE]
  for(i in seq_len(k-1)){
    out<-out+tmp[,rep.int(i+1,k),,drop=FALSE]*inveigvec2[rep.int(i+1,k),,,drop=FALSE]
  }
  out
}


microbenchmark::microbenchmark(foo(1)) #much faster, but does it work?
testy2<-foo(10.01)
plot(as.vector(testy)~as.vector(testy2))
abline(0,1)
summary(lm(as.vector(testy)~as.vector(testy2)))
cor(as.vector(testy),as.vector(testy2)) #looks really good, actually!
#and cuts down time by factor greater than 10!
#I think I'm ready to make a likelihood function based on this...
microbenchmark::microbenchmark({for(i in seq_len(nt-1)){
  out1<-C*(out1%*%P)
  out2<-C*(out2%*%P)
}})
#DRAMATIC spped up by the looks of it!
#but how bad is the bottleneck?
microbenchmark::microbenchmark(apply(test,3,eigen))
#pretty bad--it will be interesting to see whether it saves time or not...
#definitely will for larger problems, I reckon
#also might be good to see if there are any patterns here
microbenchmark::microbenchmark(lapply(eigvec,solve))
plot(eigval[1,],type='l')
plot(eigvec[1,2,],type='l')
matplot(t(rbind(eigvec[,1,],eigvec[,2,])),type='l')
#definitely some sort of pattern going on, but hard to tell what
#basically seems to transition to a matrix with 1s on the off-diagonal...interesting...
#eigenvalues certaintly have a pattern
#could also cut down on redundancy since it's all symmetric...

#now given a starting point...
gauss<-gauss.DFT(nx,dx,x0=xx[1])
out<-init<-cbind(gauss(sample(xx,1),1),gauss(sample(xx,1),0.1))
matplot(abs(out),type='l')
out3<-Re(apply(init[,1]*out1+init[,2]*out2,2,fftw::IFFT))
matplot(out3,type='l')
tesy2<-Re(apply(rep(init[,1],each=2)*testy[1,,]+rep(init[,2],each=2)*testy[2,,],1,fftw::IFFT))
tesy2<-sweep(testy,c(1,3),t(init),'*')
matplot(tesy2,type='l')
test[1,,]
matplot(t(testy[1,,]),type='l')
matplot(out1,type='l')

#simple solution to aliasing would just be to round t...
testy2<-foo(1.005)
testy3<-Re(apply(testy2,c(1,2),IFFT))
testy3[testy3<0]<- -testy3[testy3<0]
matplot(testy3[,2,],type='l')
plot(testy2[1,1,],type='l') #intuitively, the plots do seem "out of phase", but in a way that's hard to describe...
#yeah, I think some kind of rounding to nearest time point would be easiest! There doesn't seem to be an easy fix--it gets really bad as
#you approach the midpoint between samples (i.e., .05 if dx is .1)

tmp.test<-abs(apply(foo(20.05)[,1,],1,fftw::IFFT))
plot(Re(foo(1)[1,1,]),type='l')
lines(Re(foo(1.05)[1,1,]),col="red")
tmp<-foo(1.05)[1,1,]
tmp[1025:2048]<- -1*tmp[1025:2048]
lines(Re(tmp),col="blue")
plot(Re(foo(1.1)[1,1,]),type='l')
plot(Re(foo(1.15)[1,1,]),type='l')
#what happens if you just multiply right hand side by -1?

#works, but probably only for half-steps, I imagine...
#yeah, need some way of systematically adjusting right hand side...
tmp.test<-apply(foo(10.1)[,1,],1,fftw::IFFT)
matplot(abs(tmp.test),type='l')
tmp.test2<-foo(20.05)
tmp.test2[,,1025:2048]<- -1*tmp.test2[,,1025:2048]
tmp.test2<-apply(tmp.test2[,1,],1,fftw::IFFT)
matplot(abs(tmp.test2),type="l")
tmp.test2<-foo(20.05)
tmp.test2[,,1025:2048]<- -1*tmp.test2[,,1025:2048]
tmp.test2<-apply(tmp.test2[,1,],1,fftw::IFFT)

plot(foo(1)[1,1,])
plot(foo(1.05)[1,1,])
plot(foo(1.1)[1,1,])
plot(foo(1.15)[1,1,])
plot(foo(1.2)[1,1,])
tmp<-foo(1.05)[1,1,]
tmp[1025:2048]<- -1*tmp[1025:2048]
plot(tmp,
     col=rep(1:2,each=1025))
plot(foo(1.05)[1,1,],
     col=rep(1:2,each=1025))
plot(Re(tmp))
plot(Re(foo(1.05)[1,1,]))
#ah! I get it--you basically have to rotate the right-hand portion by a certain amount
plot(foo(10.05)[2,2,],
     col=rep(1:2,each=1025))
#so, a bit more complicated because the whole thing is kinda "squashed" rotationally
#seems like only the diagonal bits are affected by this for now
#(may be a parameter space thing, though?)
#need to find the amount to rotate by then stretch graph rotationally...(I think)
tmp<-foo(0.02)[,,c(1,2,2048)]
plot(foo(10.02)[1,1,],type='l')
tmp1<-Arg(tmp[1,1,2])
tmp2<-Arg(tmp[1,1,3])
if(tmp1>0){
  if(tmp2>0){
    tmp3<-2*pi-tmp2-abs(tmp1)
  }else{
    tmp3<- -tmp2-abs(tmp1)
  }
}else{
  if(tmp2>0){
    tmp3<- -tmp2+abs(tmp1)
  }else{
    tmp3<- -2*pi-tmp2+abs(tmp1)
  }
}
cor.fac<-exp(seq(0,tmp3,length.out=2048)*1i)
#well here goes nothing...
plot(foo(10.02)[1,1,]*cor.fac,col=rep(1:2,each=1024),type='l')
#seems to approximately work...
plot(abs(apply(foo(10)[,1,],1,fftw::IFFT))[,1],type="l")
lines(abs(apply(foo(10.02)[,1,],1,fftw::IFFT)*rep(cor.fac,each=2))[,1],col="red")
lines(abs(apply(foo(10.1)[,1,],1,fftw::IFFT))[,1],type="l",col="blue")
#actually seems to be pretty much exact (not always--see below)!
#pretty intuitive:
# - take difference between angle of last point in DFT and the complex conjugate of 2nd point
# - (idea being that graph of DFT should wrap around origin an integer number of times)
# - stretch graph rotationally to wrap around integer number of times...

#Gets a bit weird over long timescales...
#I wonder if it has to do with radial stretch as well? Might as well try...Nope, that wasn't it!
#Unfortunately, it might have to do with HOW MANY TIMES the graph is supposed to wrap around central point...
#Which, over long timescales, can't really be assessed because of underflow
#Holy smokes, correction factors work for the same offsets modulo t*mu/dx, I think...
#Probably just because the graph starts to go around in circles...
#But how to calculate that?
#Cheap way is to just focus on 0+modulus...
#But maybe...
#Yeah, it seems to be that deviations in the argument of the 2nd entry get wonky over long time intervals
#Need to find exact formula...
#Oh, I get it! It's because it matters which way you need to rotate
#So an exact formula is probably possible, but it's tricky because it will be hard to figure out the phase in arbitray cases where
##drifts differ between states and aren't themselves integers...
#I came up with a solution that kinda works...
#Ah, I think I get it--at a certain point, graph winds around so rapidly that you can't be sure how many revolutions ahead the next timepoint is...
#Not really sure how to prevent that, though...
#That isn't it, but I do think it has something to do with the fact that graph wraps around multiple times

#All this has lead to the insight that something odd is happening with the off-diagonal entries--they have no imaginary component for some reason...even though they really should, I think
eigvec[,,3]%*%diag(exp(10*eigval[,3]))%*%inveigvec[,,3]
#Imaginary components cancel each other out? But why?
#Maybe just not interpreting this correctly?
#Wait, no, that actually makes total sense with equal transition rates and equal drifts
#^Yeah, just a quirk of the parameters

#So after some experimenting, shockingly, aliasing only seems to occur for diagonal transition densities and is based on the drift
##parameter for that diagonal entry...
#So, the correction factors do seem able to be computed analytically, but I have no idea why
#Ah, alright--so it doesn't seem to work the more arbitrary the numbers get :(
#In the test below, neither the "brute force" nor the "analytical" methods work

tmp<-dx/c(1,-0.5)
t<-40.095
plot(foo(t)[1,1,],type='l')
plot(foo(t)[2,2,],type='l')
tmp2<-t%%tmp/abs(tmp)
tmp2[tmp2>0]<-1-tmp2[tmp2>0]-1/nx
tmp2[tmp2<0]<-tmp2[tmp2<0]+1/nx
tmp2<- -tmp2
cor.fac<-lapply(tmp2,function(ii) exp(seq(0,ii*2*pi,length.out=2048)*1i))
plot(foo(t)[1,1,]*cor.fac[[1]])
plot(foo(t)[2,2,]*cor.fac[[2]])
plot(abs(apply(foo(40)[,1,],1,fftw::IFFT))[,1],type="l")
lines(abs(apply(foo(t)[,1,]*rep(cor.fac[[1]],each=2),1,fftw::IFFT))[,1],col="red")
lines(abs(apply(foo(40.1)[,1,],1,fftw::IFFT))[,1],type="l",col="blue")

plot(abs(apply(foo(40)[,2,],1,fftw::IFFT))[,2],type="l")
lines(abs(apply(foo(t)[,2,]*rep(cor.fac[[2]],each=2),1,fftw::IFFT))[,2],col="red")
lines(abs(apply(foo(40.2)[,2,],1,fftw::IFFT))[,2],type="l",col="blue")

#Damn, there must be a way, but I'm not sure what it is...
#Maybe not--it worked for the simple case, but the more complex case just seems impossible...
#Might be a weird underflow thing--end points not quite aligned...but then again, they don't have to be...
#Seems to all come back to things wrapping around multiple times
#But even when it doesn't it fails...
#OH SHIT WAIT NO--it works!!!
#Amazing, so now you can correct the characteristic functions and infer trends despite sampling issues!!!
#It all really does boil down to making sure the characteristic function completes an integer number of revolutions
#I think you're basically just rounding things to be centered on nearest grid point
# - an improvement then might be to round the number of revolutions "down" when it's closer

#this algorithm seems to work for rounding to nearest integer revolution
#starts to fail for low values of t--still trying to work out solution
#it seems like it might be best wrap around circle by default...
#there are some really tricky edge cases to handle otherwise
# tmp<-dx/c(1,-0.5)
# t<-2.05
# plot(foo(t)[1,1,],col=rep(1:2,each=1024))
# plot(foo(t)[2,2,],col=rep(1:2,each=1024))
# tmp2<-t%%tmp/tmp
# tmp3<- tmp2-1
# inds<-abs(tmp3)<abs(tmp2)
# tmp2[inds]<-tmp3[inds]
# cor.fac<-lapply(tmp2,function(ii) exp(seq(0,ii*2*pi,length.out=2049)*1i)[-2049])
tmp<-dx/c(1,-0.5)
t<-2.05
plot(foo(t)[1,1,],col=rep(1:2,each=1024))
plot(foo(t)[2,2,],col=rep(1:2,each=1024))
tmp2<- -(t%%tmp/abs(tmp))
tmp2[tmp>0]<- -(1+tmp2)[tmp>0]
cor.fac<-lapply(tmp2,function(ii) exp(seq(0,ii*2*pi,length.out=2049)*1i)[-2049])


plot(foo(t)[1,1,]*cor.fac[[1]])
plot(foo(t)[2,2,]*cor.fac[[2]]) #makes a cardioid shape rather than a drop shape, as it should...
#also way too fat
plot(foo(t)[2,2,],col=rep(1:2,each=1024))
plot(foo(0.2)[2,2,],col=rep(1:2,each=1024))
#seems to ultimately stem from distances from origin being too high...
plot(Arg(foo(t)[2,2,]))
plot(Arg(foo(0.2)[2,2,]))
#appears to be a uniquely small time interval problem--I wonder why...
#oh my god! I get it--remember it's the inverse of the spatial domain--the distribution's too fat!
#so the argument shouldn't completely wrap around yet? Weird...
#ooh wait!
#damn, doesn't work...
#The one below is the closest you've gotten :/ But it's not good...
test<-foo(t)[2,2,]
test[1026:2048]<-Conj(test[1024:2])
# test[1025]<-complex(real=Mod(test[1025]),0) #behaves better when this is just 0...
test[1025]<-0
plot(test)
#I wonder...
test<-foo(t)[2,2,]
test<-test*exp(-Arg(test)*1i)
#doesn't work either
test<-foo(t)[2,2,]
test[1024:2048]<- -1*test[1024:2048] #this is interesting...
#what about just...http://127.0.0.1:25331/graphics/plot_zoom_png?width=1536&height=824
test<-foo(t)[2,2,]
test[1500:2048]<-0
#one last try before http://127.0.0.1:25331/graphics/plot_zoom_png?width=1536&height=824I move on!
test<-foo(t)[2,2,]
test[1026:2048]<-test[1026:2048]*exp(-(Arg(test[1026:2048])+Arg(test[1024:2]))*1i)
test[1025]<-complex(real=Mod(test[1025]),imaginary=0)


plot(log(abs(apply(foo(2)[,2,],1,fftw::IFFT))[,1]),type="l",ylim=c(-50,0))
# lines(log(abs(apply(foo(0.1)[,2,],1,fftw::IFFT))[,2]),col="red")
lines(log(abs(apply(foo(t)[,2,],1,fftw::IFFT))[,1]),col="red")
lines(log(abs(apply(foo(t)[,2,]*rep(cor.fac[[2]],each=2),1,fftw::IFFT))[,2]),col="red")
# test2<-foo(t)[2,2,]*cor.fac[[2]]
# plot(Arg(test2))
# test2[513:1536]<-test2[513:1536]*exp(-Arg(test2[513:1536])*1i)
# lines(log(abs(fftw::IFFT(test))),col="red")
# lines(log(abs(fftw::IFFT(test2))),col="red")
lines(log(abs(apply(foo(2.2)[,2,],1,fftw::IFFT))[,1]),type="l",col="blue")

#new idea:
tmp<-dx/c(1,-0.5)
t<-0.001
plot(foo(t)[1,1,],col=rep(1:2,each=1024))
plot(foo(t)[2,2,],col=rep(1:2,each=1024))
test<-foo(t)[2,2,]
test<-test*exp(-mean(Arg(test[c(2,2048)]))*1i)
plot(test)
#doesn't really seem to work...but does give an idea of how to make lower values of t work

#yeah, so there's definitely something in trying to "symmetrize" imaginary component, but it definitely doesn't work as expected :(
#not sure what you're supposed to do...

plot(Im(foo(0.01)[2,2,]))
plot(Im(foo(0.2)[2,2,]))
plot(Im(foo(0.01)[2,2,]*cor.fac[[2]]))

#After a lot of fiddling, it might be best just to round up when t is really small...
#Seems like t must complete 1 revolution--things are fine after that
tmptmp<-Re(foo(t)[2,2,]*cor.fac[[2]])
tmptmp[tmptmp<0]<-0
tmptmp<-
test<-complex(real=tmptmp,imaginary=Im(foo(t)[2,2,]*cor.fac[[2]]))

#Yeah, rounding up/down when t is small is just about the only thing I can think of...
#I just can't see how to make this work...
#Still, figuring out what to do when t is large enough is quite an achievement!
#After exploring a bit more, I'm not sure the small t problem is practical enough to warrant worrying about all too much...
#It's annoying, but you've explored your options and really can't find another way to get the job done...
#It looks bad on the log scale, but it's practically the same on absolute scale...
#Seems to become fine once you're about a tenth of the way around the circle...though it seems highly parameter-dependent
#Still, that lends confidence to the correction not being an issue in practice

#In the end, it seems an approach where you ALWAYS round up the number of revolutions works best
#There unfortunately does not seem to be an easy way to calculate the appropriate number of revolutions judging by some of
#the funky patterns you've gotten...
#And off-diagonals, for whatever reason, are still fine
