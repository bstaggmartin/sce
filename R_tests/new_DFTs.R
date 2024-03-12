#oh my god, I think I figured out the aliasing issue...
#you just needed to make the location indices symmetric about 0, rather than relying on periodicity...
#dirac deltas (any narrow distribution really) still run into issues if they are IFFT'd before getting "spread out"

.get.base<-function(nx,dx){
  tmp<-seq(0,-pi/dx,length.out=nx+1)
  base<-1i*c(tmp,-tmp[nx:2])
  tmp<-tmp^2
  list(base,c(tmp,tmp[nx:2]))
}

#Brownian motion process
BM.DFT<-function(nx,dx){
  bases<-.get.base(nx,dx)
  out<-function(t,mu,sig2){
    exp(t*(-sig2*bases[[2]]/2+mu*bases[[1]]))
  }
}

#Compound Poisson process (normal jumps)
JN.DFT<-function(nx,dx){
  bases<-.get.base(nx,dx)
  out<-function(t,lambda,tau2,delta){
    exp(t*lambda*(exp(-tau2*bases[[2]]/2+delta*bases[[1]])-1))
  }
  out
}

#normal inverse gaussian process
#seems to work--just a complicated dependence structure between tau2 and delta
#delta^2 must be less than or equal to 1/tau2
#also lambda and tau2 don't completely correspond to rate and scale of jumps...
#turns out 1/tau2 determines the steepness of jump tails (high tau2 leads to light, normalish tails)
#and lambda determines just scales the distribution width (hence why it kinda seems like a jump rate)
#although Landis interpreted lambda and tau2 as rate and jump params, so it's probably fine...
# NIG.DFT<-function(nx,dx){
#   bases<-.get.base(nx,dx)
#   out<-function(t,lambda,tau2,delta){
#     exp(t*lambda*(sqrt(1/tau2-delta^2)-sqrt(1/tau2-(delta+bases[[1]])^2)))
#   }
#   out
# }
#let's try reparameterizing delta as a proportion between -1 and 1...
#works!
#some weirdness in that the skew parameter's influence seems to become MORE severe with smaller values of tau2...
# NIG.DFT<-function(nx,dx){
#   bases<-.get.base(nx,dx)
#   out<-function(t,lambda,tau2,delta){
#     # inv.tau2<-1/tau2
#     # inv.tau<-sqrt(inv.tau2)
#     # exp(t*lambda*(inv.tau*sqrt(1-delta^2)-sqrt(inv.tau2-(inv.tau*delta+bases[[1]])^2)))
#
#     #some simplification...
#     #inv.tau2-(inv.tau*delta+bases[[1]])^2
#     #inv.tau2-(inv.tau2*delta^2+2*inv.tau*delta*bases[[1]]+bases[[2]])
#     #inv.tau*(inv.tau-inv.tau*delta^2-2*delta*bases[[1]])-bases[[2]]
#     #inv.tau*(inv.tau*(1-delta^2)-2*delta*bases[[1]])-bases[[2]]
#
#     # inv.tau<-sqrt(1/tau2)
#     # mdelta2<-1-delta^2
#     # exp(t*lambda*(inv.tau*sqrt(mdelta2)-sqrt(inv.tau*(inv.tau*mdelta2-2*delta*bases[[1]])-bases[[2]])))
#
#     #this works about as fast as anything else and avoids having to re-square bases[[1]]
#     tau<-sqrt(tau2)
#     inv.tau<-1/tau
#     # exp(t*lambda*(inv.tau*sqrt(1-delta^2)-inv.tau*sqrt(1-(delta+tau*bases[[1]])^2)))
#     mdelta2<-1-delta^2
#     exp(t*lambda*(inv.tau*sqrt(mdelta2)-inv.tau*sqrt(mdelta2-2*delta*tau*bases[[1]]+tau2*bases[[2]])))
#     #alternative parameterization where delta is squared:
#     # abs.delta<-abs(delta)
#     # exp(t*lambda*(inv.tau*sqrt(1-abs.delta)-inv.tau*sqrt(1-abs.delta-2*sign(delta)*sqrt(abs.delta)*tau*bases[[1]]+tau2*bases[[2]])))
#   }
#   out
# }
#final function
# NIG.DFT<-function(nx,dx){
#   bases<-.get.base(nx,dx)
#   out<-function(t,lambda,tau2,delta){
#     tau<-sqrt(tau2)
#     inv.tau<-1/tau
#     mdelta2<-1-delta^2
#     exp(t*lambda*(inv.tau*sqrt(mdelta2)-inv.tau*sqrt(mdelta2-2*delta*tau*bases[[1]]+tau2*bases[[2]])))
#   }
#   out
# }
#I decided to switch lambda and tau2 ultimately...better reflect behavior
NIG.DFT<-function(nx,dx){
  bases<-.get.base(nx,dx)
  out<-function(t,lambda,tau2,delta){
    inv.lambda<-1/lambda
    mdelta2<-1-delta^2
    exp(t*tau2*(inv.lambda*sqrt(mdelta2)-inv.lambda*sqrt(mdelta2-2*delta*lambda*bases[[1]]+lambda^2*bases[[2]])))
  }
  out
}

#variance gamma process
VG.DFT<-function(nx,dx){
  bases<-.get.base(nx,dx)
  out<-function(t,lambda,tau2,delta){
    (1-delta*bases[[1]]/lambda+tau2*bases[[2]]/(2*lambda))^(-t*lambda)
  }
  out
}

#dirac delta function
DI.DFT<-function(nx,dx,x0=0){
  bases<-.get.base(nx,dx)
  out<-function(x){
    exp((x-x0)*bases[[1]])
  }
}

#normal distribution
NO.DFT<-function(nx,dx,x0=0){
  bases<-.get.base(nx,dx)
  out<-function(x,epsilon2){
    exp((x-x0)*bases[[1]]-epsilon2*bases[[2]]/2)
  }
}

nx<-1024
dx<-0.05
BM<-BM.DFT(nx,dx)
xx<-c(seq(0,dx*nx,dx),seq(-dx*(nx-1),-dx,dx))
perm.inds<-c((nx+2):(2*nx),1:(nx+1))
plot(abs(fftw::IFFT(BM(10.025,-1,1)))[perm.inds]~xx[perm.inds],type="l")
abline(v=-0.05*10.025)
plot(BM(10.025,-0.05,1)) #sweet--yeah, this makes so much sense
#it was that issue where you wanted things to symmetrically spread out!

NIG<-NIG.DFT(nx,dx)
plot((Re(fftw::IFFT(NIG(1,10,64,0.5))))[perm.inds]~xx[perm.inds],type="l")

####RANDOM OU TEST####
#holy crap it actually works!

#taking derivatives...
base<-.get.base(1024,0.1)[[1]]
plot(Re(fftw::IFFT(BM(10.025,0,1)))[perm.inds]~xx[perm.inds],type="l")
lines(Re(fftw::IFFT(BM(10.025,-1,1)*(-base)))[perm.inds]~xx[perm.inds],type="l")
lines(Re(fftw::IFFT(BM(10.025,-1,1)*(-base)^2))[perm.inds]~xx[perm.inds])
lines(Re(fftw::IFFT(BM(10.025,-1,1)*(-base)^3))[perm.inds]~xx[perm.inds])
#suggests a simple pseudospectral recipe, right?
#not as much as I thought...but I found a recipe that works
#not sure about numerical stability, though, and doesn't preserve integral of function
#figured out potentially a way to avoid repeatedly fft'ing, but still working on it
#seems highly numerically unstable and can't figure out how to shift optimum to non-zero value...
# pot<-dnorm(xx,15,3)[-((nx/2+1):(3*nx/2))]
# plot(pot[perm.inds]~xx[perm.inds])
# pot<-pot/sum(pot)
nsteps<-1000
t<-5
dt<-t/nsteps
DI<-DI.DFT(nx,dx)
T1<-BM(1,15,0.1)
dbase<-diff(base[1:2])
#Okay, so I think I figured it out
#d/dx ((x-u)P) = d/dx (xP) - d/dx (uP)
#which is the derivative of the Phat (fourier transform of P) multiplied by base minus u times P multiplied by base
#shit, numerically unstable, but it works...
#there must be a way to determine the optimal timestep...
#so you get these little "oscillation bubbles forming, which I feel like MUST have some kind of solution to prevent...
#I wonder if you can just "smooth out" the frequency spectrum...
#It works REALLY well though if you just find the right timestep to use...
for(i in 2:nsteps){
  T2<-BM(dt,0,9)*T1
  # tmp<-Re(fftw::IFFT(T2))
  # runs<-rle(tmp>0)
  # inds<-runs$values&(runs$lengths>1)
  # runs$values[inds]<-FALSE
  # runs$values[!inds]<-TRUE
  # tmp[inverse.rle(runs)]<-1e-16
  # tmp[tmp<1e-16]<-1e-16
  # tmp<-tmp[-((nx/2+1):(3*nx/2))]

  #so, in the frequency domain, basically involves scaling the function by a bit
  #and convolving the derivative of the function with...well, a line?
  #is there a closed formula for that?
  #ah, the problem is that a line isn't necessarily a line in the frequency domain...
  #does have a formula, not sure if it's worth it yet...
  #so it becomes a convolution between the first derivative of dirac delta and the fourier transform of dtmp...
  #which, I think...should theoretically become the second derivative of the fourier transform of dtmp???


  # dtmp<-Re(fftw::IFFT(T2*-base))
  # dtmp<-dtmp[-((nx/2+1):(3*nx/2))]
  # tmp2<-tmp+dt*0.5*(tmp+(xx[-((nx/2+1):(3*nx/2))]+5)*dtmp) #???
  # tmp2<-c(tmp2[1:(nx/2)],rep(0,nx),tmp2[(nx/2+1):nx])
  # dT2<-T2*-base
  # ddT2<-c(dT2[1]-dT2[2*nx],diff(dT2))/dbase/2

  # dT2.test<-fftw::FFT(fftw::IFFT(T2)*xx)
  #little test to prevent numerical errors...
  #didn't work...
  # T2[(nx/2+1):(3*nx/2)]<-0

  #I also tried averaging the across x derivatives between the previous and last timestep and it did little to improve things...
  dT2<-(c(T2[-1],T2[1])-c(T2[2*nx],T2[-(2*nx)]))/(4*dbase)
  #forward diffs? Nah, central is better...
  # dT2<-(T2-c(T2[2*nx],T2[-(2*nx)]))/(2*dbase)
  # test<-Re(fftw::IFFT(T2)*xx)
  # plot(diff(test)/dx)
  # plot(Re(fftw::IFFT(dT2*-2*base)))
  T1<-T2+dt*0.5*((dT2+10*T2)*-2*base)
  # plot(Re(fftw::IFFT(T1))[perm.inds]~xx[perm.inds],type="l")

  # #one more test...is it more numerically stable to just use FFT?
  # #less so, actually!
  # dT2<-fftw::FFT(fftw::IFFT(T2)*xx)
  # T1<-T2+dt*0.5*((dT2+10*T2)*-2*base)
  # plot(Re(fftw::IFFT(T1))[perm.inds]~xx[perm.inds],type="l")

  # #so what about this?
  # #doesn't improve things either
  # dT2<-(c(T2[-1],T2[1])-c(T2[2*nx],T2[-(2*nx)]))/(4*dbase)
  # dT22<-fftw::IFFT(T2)
  # dT22<-fftw::FFT((c(dT22[-1],dT22[1])-c(dT22[2*nx],dT22[-(2*nx)]))/(2*dx))
  # T1<-T2+dt*0.5*(dT2*-2*base+10*dT22)
  # plot(Re(fftw::IFFT(T1))[perm.inds]~xx[perm.inds],type="l")

  # plot(Im(dT2))
  # plot(Im(fftw::FFT(fftw::IFFT(dT2)*-1i*xx))) #seems like this works...
  # plot(Im(diff(dT2))/dx) #seems like dividing by dx is the right thing to do...
  # plot(cumsum(Im(fftw::FFT(fftw::IFFT(dT2)*-1i*xx)))*dx)
  #hmmm...
  # tmp3<-Re(fftw::IFFT(tmp3))[-((nx/2+1):(3*nx/2))] #holy shit, that might actually work...
  # T1<-(1+dt*0.5)*(T2+c(dT2[1]-dT2[2*nx],diff(dT2)))
  # T1<-(1+dt*0.5)*(T2+fftw::FFT(fftw::IFFT(dT2)*-1i*xx)) #so something interesting going on here, but unsure...something's going very wrong any time you try iterating...
  # T1<-T2+dt*0.2*(T2+ddT2)
  # plot(Re(fftw::IFFT(T1))[perm.inds]~xx[perm.inds],type="l")
  #It works, but is numerically unstable...interesting concept though...

  # runs<-rle(tmp>0)
  # inds<-runs$values&(runs$lengths>1)
  # runs$values[inds]<-FALSE
  # runs$values[!inds]<-TRUE
  # tmp[inverse.rle(runs)]<-1e-16
  # tmp[tmp<1e-16]<-1e-16
  # plot(tmp[perm.inds]~xx[perm.inds],type="l")
  # T1<-fftw::FFT(tmp)
  if(i==2|i%%100==0){
    plot(Re(fftw::IFFT(T1)[perm.inds])~xx[perm.inds],type="l")
  }
}

#line fourier transform test
tmp<-rep(0,2*nx)
tmp[c(2,2*nx-1)]<-c(-1,1)/dx
test<-fftw::IFFT(tmp)
plot(Im(test),type="l")

tmp<-seq(min(xx),max(xx),length.out=2*nx)-5
test<-fftw::FFT(tmp)
plot(abs(test))

JN<-JN.DFT(nx,dx)
plot(abs(fftw::IFFT(JN(10.01,0.1,3,-1)))[perm.inds]~xx[perm.inds],type="l")
plot(JN(10.01,0.1,3,-1))

NIG<-NIG.DFT(nx,dx)
plot((abs(fftw::IFFT(NIG(1,1,16,0))))[perm.inds]~xx[perm.inds],type="l")

VG<-VG.DFT(nx,dx)
plot((Re(fftw::IFFT(VG(10,1,1,0))))[perm.inds]~xx[perm.inds],type="l")

DI<-DI.DFT(nx,dx,0)
plot(log(Re(fftw::IFFT(DI(-50.1)))[perm.inds])~xx[perm.inds])

NO<-NO.DFT(nx,dx,0)
plot(log(Re(fftw::IFFT(NO(-50,0.001)))[perm.inds])~xx[perm.inds],type="l")
tmp<-dnorm(xx,-50.1,sqrt(0.1),log=TRUE)
tmp<-tmp-log(sum(exp(tmp)))
lines(tmp[perm.inds]~xx[perm.inds])
plot(log(abs(fftw::IFFT(NO(-50.05,0.05)))[perm.inds])~xx[perm.inds],type="l",xlim=c(-60,-30))
tmp<-dnorm(xx,-50.05,sqrt(0.1),log=TRUE)
tmp<-tmp-log(sum(exp(tmp)))
lines(tmp[perm.inds]~xx[perm.inds])
#so dirac DFT is still vulnerable to aliasing issues, but normal DFT, crazily enough, does not seem to be...

plot(log(abs(fftw::IFFT(BM(11,-0.05,0.1)))[perm.inds])~xx[perm.inds],type="l")
plot(log(abs(fftw::IFFT(BM(11,-1.1,0.001)))[perm.inds])~xx[perm.inds],type="l")

#now you only get problems when the distribution is not sufficiently "spread out"
#this is exacerbated by things being off-period, but only a little bit
#Estimating low rates/jumps will probably be more of an issue in practice than location stuff now...

plot(log(abs(fftw::IFFT(NO(-50.1,0.01)))[perm.inds])~xx[perm.inds],type="l")
plot(log(abs(fftw::IFFT(BM(10,-1.23,0.1)))[perm.inds])~xx[perm.inds],type="l")
plot(log(abs(fftw::IFFT(NO(-50.05,0.001)*BM(10,-1.23,0.1)))[perm.inds])~xx[perm.inds],type="l")
plot(log(abs(fftw::IFFT(DI(-50.05)*BM(10,-1.23,0.1)))[perm.inds])~xx[perm.inds],type="l")

#Ooh! But these situations do seem to get "rescued" by convolution with a more spread out FFT
#So it will probably be okay in many situations--it's only a worry at the tips really
#May also want to look into rounding up values below exp(-40)--they seem to be weird aliasing artifacts

tmp<-NO(-50.1,0.05)
plot(Re(tmp)) #wish there was some kind of "nearest neighbors" thing I could employ here...
#some kind of filtering might work, but it doesn't seem to fix the main issues...
tmp[1025]<-0
plot(log(abs(fftw::IFFT(tmp))),type='l')
tmp[1025+(-100:100)]<-0
plot(log(abs(fftw::IFFT(tmp))),type='l')
#this doesn't work at all
# tmp[1025+(-100:99)]<-0+1i*Im(tmp[1025+(-100:99)])
# tmp[1025+(-100:99)]<-Re(tmp[1025+(-100:99)])+0i
# plot(log(abs(fftw::IFFT(tmp))),type='l')
tmp[c(1:101,1948:2048)]<-0
plot(log(abs(fftw::IFFT(tmp))),type='l')
tmp[1025+(-1000:1000)]<-0
plot(log(abs(fftw::IFFT(tmp))),type='l')
#yeah, if anything, leaving things unfiltered seems to be the best option, which sucks...
#oh my god, wait, the sinc function might work...

RE.DFT<-function(nx,dx){
  tmp<-seq(0,pi/dx,length.out=nx+1)
  base<-c(tmp,-tmp[nx:2])
  out<-function(width){
    tmp<-2*sin(width*base/2)/(width*base)
    tmp[1]<-1
    tmp
  }
}

#some ringing artifacts, but pretty good
RE<-RE.DFT(nx,dx)
plot(abs(fftw::IFFT(RE(100)))[perm.inds]~xx[perm.inds],type='l')

#Oh nope, you're supposed to multiple by rectangular function...
tmp<-NO(-50.1,0.001)
plot(log(abs(fftw::IFFT(tmp*RE(10)))),type='l')

#which kinda just gets you back to the beginning...
get.rect<-function(width,nx,dx){
  tmp<-round(width/(2*dx))
  c(rep(1/(2*tmp+1),1+tmp),rep(0,nx-2*tmp-1),rep(1/(2*tmp+1),tmp))
}
get.rect<-function(width,nx,dx){
  tmp<-round(width/(2*dx))
  c(rep(1,1+tmp),rep(0,nx-2*tmp-1),rep(1,tmp))
}
plot(get.rect(100,nx,dx),type="l")
sum(get.rect(100,nx,dx))
plot(abs(fftw::IFFT(tmp*get.rect(102,nx,dx))),type='l')
plot(abs(fftw::IFFT(tmp)),type='l')

#yeah, doesn't really seem to work well--at least, not any better than the unfiltered result

#Think I figured it out! Will potentially make algorithm annoying, but I think it's better than nothing:
#Use Re() instead of abs()--seems to naturally find cut-off points for when estimated densities are reliable

DI<-DI.DFT(nx,dx,0)
tmp<-Re(fftw::IFFT(DI(-50.1)))
plot(log(tmp)~xx)
plot(log(tmp)[perm.inds]~xx[perm.inds],type='l')
plot(tmp>0~xx) #doesn't seem like much, but they often have a single "island"
test<-rle(tmp>0)
plot(test$lengths[test$values]) #only 1 "island" of positive values!

tmp<-Re(fftw::IFFT(NO(-50.5,0.1)))
plot(log(tmp)~xx)
plot(log(tmp)[perm.inds]~xx[perm.inds],type='l')
plot(tmp>0~xx) #doesn't seem like much, but they often have a single "island"
test<-rle(tmp>0)
plot(test$lengths[test$values]) #only 1 substantial "island" of positive values in most cases!

#recipe for isolating the island:
test<-rle(tmp>0)
ind<-which(test$values)
ind<-ind[which.max(test$lengths[ind])]
test$values[ind]<-FALSE
test$values[-ind]<-TRUE
tmp[inverse.rle(test)]<-0

plot(log(tmp[perm.inds])~xx[perm.inds],xlim=c(-40,-60),type='l')
tmp2<-dnorm(xx,-50.1,sqrt(0.1),log=TRUE)
tmp2<-tmp2-log(sum(exp(tmp2)))
lines(tmp2[perm.inds]~xx[perm.inds]) #perfect!
#I guess it might run into issues with bimodal distributions, now that I think on it...

#How different when you take all islands?
test<-rle(tmp>0)
inds<-test$values&(test$lengths>1)
test$values[inds]<-FALSE
test$values[!inds]<-TRUE
tmp[inverse.rle(test)]<-1e-16
tmp[tmp<1e-16]<-1e-16

#overall seems to be a suitably minimal effect
plot(log(tmp))
plot(log(tmp[perm.inds])~xx[perm.inds],type='l')
tmp2<-dnorm(xx,-50.5,sqrt(2),log=TRUE)
tmp2<-tmp2-log(sum(exp(tmp2)))
lines(tmp2[perm.inds]~xx[perm.inds])
# lines(log(abs(fftw::IFFT(NO(-50.1,0.1))))[perm.inds]~xx[perm.inds])
points(log(Re(fftw::IFFT(NO(-50.5,2))))[perm.inds]~xx[perm.inds])
#when you think about all the other stuff going on in these algorithms, I feel like the chances of the
#"false positive" islands screwing with results is quite minimal
#I think there's too much weird variation in the tail behavior otherwise--this is the "safe" way to go
#(although it certain makes me worry about how frequent "holes" in the probability surface will be...)
#(one work-around might be to replace 0s with suitably low numbers--like 1e-16)
#(I would just implement the cleaning algorithm both ways and see which one works better...)
#Yeah, looking into it, a two-step process would be best:
# - Round any isolated values with positive real part
# - Round any values with real part below 1e-16 (~-37 on log scale)
#These can either be rounded to 1e-16 (incorrect tail behavior, but prevent catastrophic underflow)...
#...or 0 (more in line with tail behavior, but underflow potentially catastrophic)
#basically, you're looking to keep any values above 1e-16 with at least one neighbor above 1e-16
#that can be done quite easily...
tmp<-Re(fftw::IFFT(DI(-100.49)))
ll<-c(2*nx,1:(2*nx-1))
rr<-c(2:(2*nx),1)
inds<-tmp<1e-16
# inds<-inds|(inds[ll]&inds[rr]&tmp<0.5*max(tmp))
#replace underflowed values AND low-ish values sorrounded by underflowed values with 1e-16 (or maybe 0)
tmp[inds|(inds[ll]&inds[rr]&tmp<0.5*max(tmp))]<-1e-16
#rescale to make sure distribution still sums to unity
# inds<-!inds
# tmp[inds]<-tmp[inds]/sum(tmp[inds])
tmp<-tmp/sum(tmp)
#do you allow the rescaling step to rescale the 1e-16?
#yeah, I think so--it barely changes anything at that point...
plot(log(tmp))
plot(tmp)
#a bit faster a more conceptually straight-forward
#only issue is eliminating true islands of length 1 (especially clean dirac functions, etc.)
#just adding an extra condition to keep anything greater than half the max height seems appropriate...
#seems to work, and decently fast too!
#I think this is the way
#Oh, but should rescale to sum to 1 too...



library(inline)
library(Rcpp)
library(RcppArmadillo)
Rcpp::cppFunction("
arma::cx_cube vec_exp(arma::cx_cube x) {

  arma::cx_cube res = x;

  res.each_slice( [](arma::cx_mat& X) { X = arma::expmat(X); } );

  return res;
}
",
depends="RcppArmadillo")
Rcpp::cppFunction("
List vec_eig(arma::cx_cube x) {

  arma::cx_cube eigvecs(x.n_rows, x.n_cols, x.n_slices);
  arma::cx_cube inveigvecs = eigvecs;
  arma::cx_mat eigvals(eigvecs.n_rows, eigvecs.n_slices);
  arma::cx_mat eigvec(x.n_rows, x.n_cols);
  arma::cx_vec eigval(x.n_rows);

  for(int i=0; i<x.n_slices; i++) {
    eig_gen(eigval, eigvec, x.slice(i));
    eigvals.col(i) = eigval;
    eigvecs.slice(i) = eigvec;

  }

  return List::create(eigvals, eigvecs);
}
",
depends="RcppArmadillo")

BM.DFT2<-function(nx,dx){
  bases<-.get.base(nx,dx)
  out<-function(mu,sig2){
    -sig2*bases[[2]]/2+mu*bases[[1]]
  }
}

nx<-1024
dx<-0.1

test<-array(c(-2,1,2,-1),c(2,2,2*nx))
tmp<-BM.DFT2(nx,dx)
test[1,1,]<-test[1,1,]+tmp(2,10)
test[2,2,]<-test[2,2,]+tmp(-4,100)

eig<-lapply(asplit(test,3),eigen)
vecs<-lapply(eig,'[[',"vectors")
invvecs<-lapply(vecs,solve)
vecs<-array(unlist(vecs,use.names=FALSE),c(2,2,2*nx))
invvecs<-array(unlist(invvecs,use.names=FALSE),c(2,2,2*nx))
vals<-do.call(cbind,lapply(eig,'[[',"values"))
foo<-function(t){
  tmp<-sweep(vecs,c(2,3),exp(t*vals),'*',check.margin=FALSE)
  out<-tmp[,rep.int(1,2),,drop=FALSE]*invvecs[rep.int(1,2),,,drop=FALSE]
  for(i in seq_len(2-1)){
    out<-out+tmp[,rep.int(i+1,2),,drop=FALSE]*invvecs[rep.int(i+1,2),,,drop=FALSE]
  }
  out
}
foo2<-function(t,x,k){
  tmp<-invvecs[,1,]*x[rep.int(1,k),]
  for(i in seq_len(k-1)){
    tmp<-tmp+invvecs[,i+1,]*x[rep.int(i+1,k),]
  }
  tmp<-exp(t*vals)*tmp
  out<-vecs[,1,]*tmp[rep.int(1,k),]
  for(i in seq_len(k-1)){
    out<-out+vecs[,i+1,]*tmp[rep.int(i+1,k),]
  }
  out
}

#still way faster to do eigen decomp...
#but you could always use eigen decomp with c?
microbenchmark::microbenchmark(foo(1))
rcpp_res<-vec_solve(test)
r_res<-foo(1)
plot(as.vector(log(abs(r_res)))~as.vector(log(abs(rcpp_res)))) #may get better results with rcpp...
plot(log(abs(Re(rcpp_res))),type="l")
plot(log(abs(Re(r_res))),type="l")
test.ifft<-fftw::IFFT(0.5*rcpp_res[1,1,]+0.5*rcpp_res[1,2,])
plot(log(Re(test.ifft)),type="l")
test.ifft<-fftw::IFFT(0.5*r_res[1,1,]+0.5*r_res[1,2,])
plot(log(Re(test.ifft)),type="l")
#interesting! If anything, actually get better results with eigendecomp--wonder why...
Rcpp::cppFunction('
List vec_eig(arma::cx_cube x) {

  arma::cx_mat eigvals(x.n_rows, x.n_slices);
  arma::cx_vec eigval(x.n_rows);
  arma::cx_cube eigvecs(x.n_rows, x.n_cols, x.n_slices);
  arma::cx_cube inveigvecs = eigvecs;
  arma::cx_mat eigvec(x.n_rows, x.n_cols);


  for(int i=0; i<x.n_slices; i++) {
    eig_gen(eigval, eigvec, x.slice(i));
    eigvals.col(i) = eigval;
    eigvecs.slice(i) = eigvec;
    inveigvecs.slice(i) = arma::inv(eigvec);
  }

  return List::create(eigvals, eigvecs, inveigvecs);
}
',
depends="RcppArmadillo")

#naive version--not that fast
# Rcpp::cppFunction("
# arma::cx_cube vec_exp(double t, arma::cx_mat eigvals, arma::cx_cube eigvecs, arma::cx_cube inveigvecs) {
#
#   eigvals = exp(t * eigvals);
#
#   for(int i=0; i<eigvecs.n_slices; i++) {
#     eigvecs.slice(i) = eigvecs.slice(i) * diagmat(eigvals.col(i)) * inveigvecs.slice(i);
#   }
#
#   return eigvecs;
# }
# ",
# depends="RcppArmadillo")

Rcpp::cppFunction("
arma::cx_cube vec_exp(double t, arma::cx_mat eigvals, arma::cx_cube eigvecs, arma::cx_cube inveigvecs) {

  arma::cx_mat vals(eigvals.n_rows, eigvals.n_cols);
  arma::cx_cube vecs = eigvecs;

  vals = exp(t * eigvals);
  for(int i=0; i<eigvals.n_rows; i++) {
    vecs.row(i) = vecs.row_as_mat(i).t() % vals;
  }

  arma::cx_cube out(eigvecs.n_rows, eigvecs.n_cols, eigvecs.n_slices);
  arma::cx_cube tmp1 = out;
  arma::cx_cube tmp2 = out;
  for(int i=0; i<eigvals.n_rows; i++) {
    for(int j=0; j<eigvals.n_rows; j++) {
      tmp1.col(j) = vecs.col(i);
      tmp2.row(j) = inveigvecs.row(i);
    }
    out = out + tmp1 % tmp2;
  }

  return out;
}
",
depends="RcppArmadillo")

microbenchmark::microbenchmark(lapply(asplit(test,3),eigen),times=10)
microbenchmark::microbenchmark(vec_eig(test),times=10) #soooo much faster!!!
#honestly, you may just want to implement the entire likelihood function in C...
#once you calculate the Q array, it's just a matter of eigendecomp, multiplication, and ffts...
tmp<-vec_eig(test)
microbenchmark::microbenchmark(foo(1))
microbenchmark::microbenchmark(vec_exp(1, tmp[[1]], tmp[[2]], tmp[[3]])) #a bit faster...
#ooh, can speed up even faster by pre-multiplying inverse cube with weight matrix
vec_exp(1, tmp[[1]], tmp[[2]], tmp[[3]])[,,500]
foo(1)[,,500]
tmp[[2]][,,500]
vecs[,,500]
tmp[[1]][,500]
vals[,500]

sweep(tmp,2,diag(tmp2),"*")%*%rowSums(sweep(tmp3,2,tmp4,"*"))

#this is much faster, but seems numerically unstable...
#only because of matrix exponential, as it turns out!
#eigendecomposition is surprisingly fine
#even manually calculating complex exponential results in problems...

Rcpp::cppFunction("
arma::cx_mat vec_exp_mult(double t, arma::cx_mat x, arma::cx_mat eigvals, arma::cx_cube eigvecs, arma::cx_cube inveigvecs) {

  arma::cx_mat out1(eigvals.n_rows, eigvals.n_cols);
  arma::cx_mat tmp_exp(eigvals.n_rows, eigvals.n_cols);
  arma::cx_mat out2(eigvals.n_rows, eigvals.n_cols);

  for(int i=0; i<eigvals.n_rows; i++) {
    out1.row(i) = arma::sum(inveigvecs.row_as_mat(i).t() % x);
  }
  // out1 %= arma::exp(t * eigvals);
  tmp_exp.set_real(exp(t * real(eigvals)) % cos(t * imag(eigvals)));
  tmp_exp.set_imag(exp(t * real(eigvals)) % sin(t * imag(eigvals)));
  out1 %= tmp_exp;
  for(int i=0; i<eigvals.n_rows; i++) {
    out2.row(i) = arma::sum(eigvecs.row_as_mat(i).t() % out1);
  }

  return out2;
}
",
depends="RcppArmadillo")

tmp<-vec_eig(test)
xx<-fftw::FFT(c(rep(0,nx),1,rep(0,nx-1)))
xx<-rbind(xx,xx)
res<-vec_exp_mult(1, xx, vals, vecs, invvecs)
plot((Re(fftw::IFFT(res[1,]))),type="l")
plot((Re(fftw::IFFT(res[2,]))),type="l")
microbenchmark::microbenchmark(vec_exp_mult(1, xx, tmp[[1]], tmp[[2]], tmp[[3]]))

tmpy<-foo(1)
tmpy<-do.call(cbind,lapply(seq_len(2*nx),function(ii) tmpy[,,ii]%*%xx[,ii]))
plot(Re(fftw::IFFT(tmpy[1,])),type="l")
plot(Re(fftw::IFFT(tmpy[2,])),type="l")

tmpy<-vec_exp(1, tmp[[1]], tmp[[2]], tmp[[3]])
tmpy<-do.call(cbind,lapply(seq_len(2*nx),function(ii) tmpy[,,ii]%*%xx[,ii]))
plot(Re(fftw::IFFT(tmpy[1,])),type="l")
plot(Re(fftw::IFFT(tmpy[2,])),type="l") #hmmm...

tmpy<-vec_exp(1, vals, vecs, invvecs)
tmpy<-do.call(cbind,lapply(seq_len(2*nx),function(ii) tmpy[,,ii]%*%xx[,ii]))
plot((Re(fftw::IFFT(tmpy[1,]))),type="l")
plot((Re(fftw::IFFT(tmpy[2,]))),type="l")

vals<-tmp[[1]]
vecs<-tmp[[2]]
invvecs<-tmp[[3]]
tmpy<-foo2(1,xx,2)

microbenchmark::microbenchmark(foo(1,xx,2))


#so R has wayyy more stable eigendecomposition
?eigen
#looks like you will need to use RcppEigen to do stabler

Rcpp::cppFunction("
List vec_eig(arma::cx_cube x) {

  arma::cx_mat eigvals(x.n_rows, x.n_slices);
  arma::cx_vec eigval(x.n_rows);
  arma::cx_cube eigvecs(x.n_rows, x.n_cols, x.n_slices);
  arma::cx_cube inveigvecs = eigvecs;
  arma::cx_mat eigvec(x.n_rows, x.n_cols);

  arma::cx_mat tmp_mat1(x.n_rows, x.n_cols);
  Eigen::MatrixXcf tmp_mat2(x.n_rows, x.n_cols);

  for(int i=0; i<x.n_slices; i++) {
    tmp_mat = x.slice(i);

    eig_gen(eigval, eigvec, x.slice(i));
    eigvals.col(i) = eigval;
    eigvecs.slice(i) = eigvec;
    inveigvecs.slice(i) = arma::inv(eigvec);
  }

  return List::create(eigvals, eigvecs, inveigvecs);
}
",
depends="RcppArmadillo")

#R offers way better numerical stability
#But you can use the right-to-left matrix multiplication approach you found to slightly speed things up...
#Shockingly, seems to have more to  do with complex multiplication than eigendecomposition...
#Suggesting that you may actually be best served by using RcppArmadillo for eigendecomposition, followed by R for multiplication
#Still saving a shitload of time--woohoo!

#Maybe something in vec_exp_mult isn't behaving as intended? It's hard to believe it would be this bad with complex numbers...
#Either that or the complex exponential function is off...
#Oh weird, it's because the complex exponential function has poor numerical precision!
#Welp, that sucks, but it doesn't seem to save much time to use armadillo on the multiplication steps anyways...something to explore later
