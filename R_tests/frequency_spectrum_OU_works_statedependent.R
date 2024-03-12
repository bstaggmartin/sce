#1/25: over the timescales investigated here, with realistic empirical parameters, the risk of boundary conditions really screwing with things seems pretty small...
#Not non-existent, but I think there's a good chance the approximations will just wash out in the final scheme!

#logistic potential?
# plot(Re(1/sinh(base*nx*dx^2*1i/pi))) #something like that
#convolution with this would approximately just look like derivative though...
#you'd have to calculate the full convolution probably for it to matter
#yeah, I think you found all the ones that make sense...otherwise you'd end up with a very dense matrix that would be difficult to solve

#one more thing--periodic linear potential...
#sum of heaveside step and linear function...right?
unit.step<-c(0.5,rep(1,(2*nx-1)/2),0.5,rep(0,(2*nx-1)/2))
step2<-c(rep(0,nx/2),0.5,rep(1,nx-1),0.5,rep(0,nx/2-1))
plot((xx+nx*dx*step2)) #hmmm...no, not quite
#you'd need to flip the sign of the step function depending on the side...
step3<-c(rep(0,nx/2),-0.5,rep(-1,nx/2-1),0,rep(1,nx/2-1),0.5,rep(0,nx/2-1))
#hmmm...easier to think of as two rect functions at that point
test<-c(rep(0,nx/2+1),rep(-1,nx/2),rep(1,nx/2),rep(0,nx/2-1))
plot(test[perm.inds])
plot(xx+nx*dx*test)
plot(Im(IFFT(xx+nx*dx*test)))
plot(Im(IFFT(xx)))
#so very similar to first derivative of dirac, but has these alternating 0s and imaginary components...
plot(Im(FFT(xx+nx*dx*test)),type="l",xlim=c(1900,2048))
plot(Im(FFT(xx)),type="l",xlim=c(1900,2048))
#Ohhh, I get it--first dirac derivative oscillates with consistent period of 2*1, while xx+test starts exhibit sinc-type stuff with the periodicity of oscillations increasing!
plot(abs(FFT(test)),type="l",xlim=c(0,100))
abline(v=seq(1,100,2))
#Oh wait, actually has consistent period of 2*2 for xx+test...
#Real component's odd though...consistent 0,1,-2 pattern
#Convolution with this would mean an odd kind of "double derivative", right?
#I feel like taking the central difference in the pentadiagonal way might be how to do it?
plot(Im(FFT((xx+nx*dx*test)*dnorm(xx,45,100))),type="l")
tmp<-FFT(dnorm(xx,45,100))
indsrr<-c(3:(2*nx),1:2)
indsll<-c(2*nx-1,2*nx,1:(2*nx-2))
plot(Im((tmp[indsrr]-tmp[indsll])/(4*dbase)),type="l")
plot(Im(FFT((xx+nx*dx*test)*dnorm(xx,45,1)))~Im((tmp[indsrr]-tmp[indsll])/(4*dbase)))
#Yeah, that does the trick it seems...
#oh never mind, not quite! Messes up when probability mass is in the padded region, it seems
#which means that it probably is a pretty complicated convolution
#possible to find some approx solution perhaps...but not obvious at this juncture
#yeah, just an absolute mess--may be worth writing things out at some point to see if you can find a helpful identity somewhere, but I'm not gonna hold my breath
#It feels like at this rate it would be easier to IFFT, enforce boundary conditions, and FFT again
#It's just not an easy convolution


plot(((xx+nx*dx*test)*dnorm(xx,45,100))[perm.inds]~xx[perm.inds],type="l")
lines(Im(IFFT((tmp[indsrr]+tmp[indsll]-tmp)/(4*dbase)))[perm.inds]~xx[perm.inds],type="l")

#What about when the function is shifted?
#well, in light of the fact that the imaginary part is complicated...who knows
#but it seems like there's some established patterns that might be useful
plot(2*(xx+nx*dx*test-5*(1+test))[perm.inds]~xx[perm.inds],type="l")
abline(v=5)
abline(h=0)
#almost, but not quite...unfortunately seems to require more complex combo of functions...
test1<-c(rep(0,nx/2+1),rep(1,nx/2),rep(0,nx-1))
test2<-c(rep(0,nx+1),rep(1,nx/2),rep(0,nx/2-1))
plot(2*(xx+nx*dx*test-5+10*test1+10*test2)[perm.inds]~xx[perm.inds],type="l") #seems to work
#so you get a combo of step function AND rectangular functions...
#which seems potentially quite problematic...
#but it may be worth exploring this thread a little further given the surprising pentadiagonal-reflective boundary condition you noticed
step4<-c(rep(0,nx/2+1),rep(1,nx),rep(0,nx/2-1))
plot(Re(FFT(step4)))
plot(Im(FFT(test)))
plot(Re(FFT(xx+nx*dx*test-10+2*10*step4))) #causes some weirdness with the real part, though the imaginary part is interestingly unaffected
plot(Re(FFT(xx+nx*dx*test-10)))
plot(Re(FFT(xx+nx*dx*test-10+2*10*step4))~Re(FFT(xx+nx*dx*test-10)),xlim=c(-10,10),ylim=c(-10,10))

plot(Im(FFT(dnorm(xx,5,50)*(xx+nx*dx*test-10+2*10*step4))),type="l")
plot(Im(FFT(dnorm(xx,5,50)*(xx+nx*dx*test-10))),type="l")
plot(Im(FFT(dnorm(xx,5,50)*(xx+nx*dx*test-10+2*10*step4)))~Im(FFT(dnorm(xx,5,50)*(xx+nx*dx*test-10))))
abline(0,1) #nearly the same, but not quite...seems like things "deflect" ever so slightly
plot(Re(FFT(dnorm(xx,5,50)*(xx+nx*dx*test-10+2*10*step4))),type="l")
plot(Re(FFT(dnorm(xx,5,50)*(xx+nx*dx*test-10))),type="l")
plot(Re(FFT(dnorm(xx,5,50)*(xx+nx*dx*test-10+2*10*step4)))~Re(FFT(dnorm(xx,5,50)*(xx+nx*dx*test-10))))
abline(0,1)

#idea: a dynamic windowing function that allows change in any sufficiently large values of the initial distribution
#then, once they can be safely rounded to 0, ignores them...
#while ensuring that any values corresponding to the narrowest stationary distribution can continue to change?
#so there's something to this, but it's not perfect...
#essentially, you want the solution to get "smoother" over time
#the irregular spikes seem pretty obvious most of the time
#the issue is how do we detect and remove these spikes? Because they can just result from the initial distribution
#so you need some way of detecting "new" spikes while not freaking out about the old ones
#which actually seems shockingly hard to do...
#oh shoot, and you're going backwards, which means the variance is unbounded...
#but that's actually kinda useful because it means the central spike of the frequency domain should always be decreasing...
#yeah, there's probably a rough way you can leverage that intuition to your advantage--as soon as a wave decreases to 0 amplitude, it should remain that way unless it receives "flux" from the other states...
#which is easily calcuable during the convolution step!
#so keep track of "hot" cells capable of changing and "cold" cells locked in place (rounded to 0? or will that fuck things up?)
#cold cells can only be rendered hot again from a hot cell flowing to a cold cell via the convolution step...
#may render the foo() part slightly innaccurate (ignoring derivative info), but that may be worth the cost...
#wait, if I'm not mistaken, all this could be implemented by simply rounding values to 0...right?
#so it really does seem like the detectable signature is the non-monotonic decay of absolute values to high negative/positive freqs
#questions: would it help with the weird "anti-optimum" stuff if you do the fourier transform from 0 to 2pi rather -pi to pi?
#how would you detect and remove non-monotonic bumps?
#yeah, it's hard because non-monotonic bumps easily result from initial conditions
#the artifacts result in rapid oscillations in the frequency domain, but it remains mysterious how to detect/prevent such oscillations from forming...
#yeah, I tried a bunch of "simple" solutions involving smoothing out the finite differences step, but nothing seems to work as well you'd hoped...
#the key to understanding the oscillations is that they result from things undershooting to negative numbers, but there's no easy way to prevent this in the frequency domain...
#the only thing that really worked was over-aggressive smoothing...

#Yeah, I think unfortunately the only way you're gonna solve this is some implicit scheme
#The numerical differentiation really seems to work quite well, and it's those over/undershoots you have to work to prevent...
#The only thing that seems to consistently work decently well is smoothing and FFT'ing at each step to ensure things remain positive...
#Even then, doesn't prevent spurious oscillations, just prevents them from blowing up

#The new implicit tridiagonal matrix method works fantastically!
#Definitely gets a bit wonky going backwards with the periodic boundary conditions, but I think you can come up with some hacks to that end
#On the more analytical end of things, an interesting posibility would be to instead use a periodic cubic potential parameterized...
#such that the main part of the potential is preserved...
#i.e.
plot(alpha1*(xx^2-5*xx)[perm.inds]~xx[perm.inds])
#then you have to find a cubic term that makes this function meet at the boundaries...
tmp<-alpha1*(xx^2-5*xx)[perm.inds]
(tmp[2]-tmp[2*nx])
plot((tmp[2]-tmp[2*nx])/2*(xx^3/(nx*dx)^3))
plot((tmp[2]-tmp[2*nx])/2*(xx[perm.inds]^3/(nx*dx)^3)+tmp) #sweet
tmp2<-(tmp[2]-tmp[2*nx])/2*(xx[perm.inds]^3/(nx*dx)^3)+tmp
plot(diff(tmp2))
plot(diff(tmp)) #hmmm...doesn't make the derivative continuous though, which would probably still cause problems
#you would need a quartic term in the potential, which would require 3rd derivatives, which would probably be unstable to try numerically...
tmp<-alpha1*(xx-5)[perm.inds]
plot((tmp[2]-tmp[2*nx])/2*(xx[perm.inds]^3/(nx*dx)^3))
plot((tmp[2]-tmp[2*nx])/2*(xx[perm.inds]^3/(nx*dx)^3)+tmp~xx[perm.inds])
abline(v=5)
tmp2<-(tmp[2]-tmp[2*nx])/2*(xx[perm.inds]^3/(nx*dx)^3)+tmp
plot((tmp[2]-tmp[2*nx])/8*(xx[perm.inds]^4/(nx*dx)^3)+alpha1*xx[perm.inds]^2/2-alpha1*5*xx[perm.inds]~xx[perm.inds])
#well, I at least get why things get stuck where they get they do
#another idea would be to use a trig curve potential, but there's not really an easy way to independently control both the periodicity and slope...is there?
#maybe fix the periodicity and modulate the amplitude to get the right shape?

#two unstable states at the midpoint between the optimum and antioptimum:
plot(alpha1*(xx^2/2-5*xx)[perm.inds]~xx[perm.inds],type="l")
lines(-alpha1*nx^2*dx^2/(4*pi^2)*cos(2*pi*(xx[perm.inds]-5)/(nx*dx))+alpha1*nx^2*dx^2/(4*pi^2)~xx[perm.inds])
#so you want it to have a derivative of alpha1 about 5...
plot(diff(alpha1*(xx^2/2-5*xx)[perm.inds]),type="l")
lines(diff(-alpha1*nx^2*dx^2/(4*pi^2)*cos(2*pi*(xx[perm.inds]-5)/(nx*dx))))
plot(alpha1*(xx[perm.inds]-5)~xx[perm.inds])
lines(alpha1*nx*dx/(2*pi)*sin(2*pi*(xx[perm.inds]-5)/(nx*dx))~xx[perm.inds])
abline(v=-nx*dx/2)
abline(v=nx*dx/2)
abline(h=0)
#this definitely might help...would theoretically help impose the proper boundary conditions

#one unstable state at the anti-optimum:
plot(alpha1*(xx^2/2-25*xx)[perm.inds]~xx[perm.inds],type="l")
lines(-alpha1*nx^2*dx^2/(pi^2)*cos(pi*(xx[perm.inds]-25)/(nx*dx))+alpha1*nx^2*dx^2/(pi^2)~xx[perm.inds])
plot(diff(alpha1*(xx^2/2-5*xx)[perm.inds]),type="l")
lines(diff(-alpha1*nx^2*dx^2/(pi^2)*cos(pi*(xx[perm.inds]-5)/(nx*dx))))
plot(alpha1*(xx[perm.inds]-10)~xx[perm.inds])
lines(alpha1*nx*dx/(pi)*sin(pi*(xx[perm.inds]-10)/(nx*dx))~xx[perm.inds])
abline(v=-nx*dx/2)
abline(v=nx*dx/2)
abline(h=0)

#if I'm reading my table of fourier transform properties correctly...using a trigonometric function-based potential should make things even easier!
#ends up looking like adding two shifted copies of fourier transform to unshifted copy
#the only trick will be figuring out how to solve permuted tridiagonal matrices...which I think is possible but I'm unsure
#yeah, that will end up being the trickiest part for sure
#welp, you done good for today I think
#oh shit no! look at the formula-->I think it will always end up being tridiagonal with 2*pi/(nx*dx)
#base-2*pi/(2*pi)/(nx*dx)=base-1/(nx*dx)
#except...wait...isn't the base purely imaginary? shoot
#so if a isn't complex it doesn't make sense from a numeric perspective, right?
#there's also the matter of the constant in the sine function, though I feel like there might be a helpful identity there...
#I feel like this isn't that complicated and I'm missing something...
test<-BM(1,0,1)
plot(Re(IFFT(test)[perm.inds])~xx[perm.inds],type="l")
plot(Im(test),type="l")

plot(Re(FFT(Re(IFFT(test))*sin(2*pi*(xx-5)/(nx*dx)))),type="l") #that seems...right?
#but yeah, that does actually work
#and, funny enough, it matches quite well with what you might expect given the derivative convolution stuff...
#nice!
#what about theta?
#seems like it affects the scale and recaptiulates a flipped version of the unshifted function in the reals
#but it leaves the shape of the imaginary part unchanged!
#definitely could figure this out with enough puzzling...

plot(Re(FFT(Re(IFFT(test))*sin(2*pi*(xx-5)/(nx*dx)))),type="l")
#ugh, I can get the real or imaginary coefficients right, but getting them both right presents a challenge...
plot(Re((exp(-5i/(nx*dx))*BM(1,dx,1)-exp(5i/(nx*dx))*BM(1,dx,1))),type="l")
Im((exp(2i*pi*5/(nx*dx))*BM(1,dx,1)-exp(-2i*pi*5/(nx*dx))*BM(1,dx,1))/(-2i))/Im(FFT(Re(IFFT(test))*sin(2*pi*(xx-5)/(nx*dx))))

#let's simplify a bit...
plot(Re(FFT(sin(2*pi*(xx-5)/(nx*dx)))),type="l")
rr<-c(0,0,1,rep(0,2*nx-3))
ll<-c(rep(0,2*nx-2),1,0)
plot(Re((exp(1i*2*pi*-5/(nx*dx))*2*nx*rr-(exp(-1i*2*pi*-5/(nx*dx))*2*nx*ll))/2i),type="l") #nice! got it!

rr<-test[c(2*nx-1,2*nx,1:(2*nx-2))]
ll<-test[c(3:(2*nx),1:2)]
plot(Im(FFT(Re(IFFT(test))*sin(2*pi*(xx-5)/(nx*dx)))),type="l")
plot(Im((exp(-1i*2*pi*5/(nx*dx))*rr-exp(1i*2*pi*5/(nx*dx))*ll)/2i),type="l") #got it!
#so you actually wanna do pi for a properly tridiagonal matrix...
rr<-test[c(2*nx,1:(2*nx-1))]
ll<-test[c(2:(2*nx),1)]
plot(Re(FFT(Re(IFFT(test))*sin(pi*(xx-5)/(nx*dx)))),type="l")
plot(Re((exp(-1i*pi*5/(nx*dx))*rr-exp(1i*pi*5/(nx*dx))*ll)/2i),type="l")
#(exp(-1i*pi*theta/(nx*dx))*rr-exp(-1i*pi*theta/(nx*dx))*ll)/2i
#interestingly, this implies a diagonally dominant hermitian matrix as long as dt*alpha is less than 1!
#oh, never mind...because base doesn't multiply the -1s of the diagonal...
#so it actually might be less stable in some ways...
#But still kinda gives you more insight into the structure of the matrix

#the constant to multiply everything by, by the way, is alpha1*nx*dx/pi



library(fftw)

nx<-1024
dx<-0.05
BM<-BM.DFT(nx,dx)
xx<-c(seq(0,dx*nx,dx),seq(-dx*(nx-1),-dx,dx))
perm.inds<-c((nx+2):(2*nx),1:(nx+1))
base<-.get.base(nx,dx)
base2<-base[[2]]
base<-base[[1]]
baseish<-base*BM(1,0,dx)
dbase<-diff(base[1:2])
inds1<-c(2:(2*nx),1)
inds2<-c(2*nx,1:(2*nx-1))
k<-2
Q<-0.02*matrix(c(-1,1,1,-1),2,2)

alpha1<-0.02
alpha2<-0.04
theta1<-3.6
theta2<- -2.1
sig21<-0.04
sig22<-0.1
#can we treat these as exponentials???
#that didn't work either... :(
#higher order differences do seem to help, but not too much...
foo1<-function(x){
  alpha1*base*(theta1*x+(x[indsll]/12+2*x[inds2]/3-2*x[indsr]/3-x[indsrr]/12)/dbase)
}
# #matrix form of linear differential operator...
# mat<-matrix(0,2*nx,2*nx)
# diag(mat)<-base*theta1#-sig21/2*base2 #the sig2 part adds in the diffusion component, which might make things more stable for implicit solving...
# mat[cbind(seq_len(2*nx),inds1)]<- -base/(2*dbase)
# mat[cbind(seq_len(2*nx),inds2)]<- base/(2*dbase)
# mat<-alpha1*mat
# mat2<-mat
# diag(mat2)<-base*theta2
# mat2[cbind(seq_len(2*nx),inds1)]<- -base/(2*dbase)
# mat2[cbind(seq_len(2*nx),inds2)]<- base/(2*dbase)
# mat2<-alpha2*mat2

#quadratic potential
a1<-alpha1*base/(2*dbase)
b1<-alpha1*base*theta1
c1<- -alpha1*base/(2*dbase)
a2<-alpha2*base/(2*dbase)
b2<-alpha2*base*theta2
c2<- -alpha2*base/(2*dbase)
a1<- -a1
b1<- -b1
c1<- -c1
a2<- -a2
b2<- -b2
c2<- -c2

#trigonometric potential
#these expression could be further simplified for sure...
#seems to yield much better behavior
#would be nice to do the more highly periodic potential, but then you have to figure out solving banded matrices...
a1<- -alpha1*nx*dx*base/pi*exp(-1i*pi*theta1/(nx*dx))/2i
b1<-rep(0,2*nx)
c1<-alpha1*nx*dx*base/pi*exp(1i*pi*theta1/(nx*dx))/2i
a2<- -alpha2*nx*dx*base/pi*exp(-1i*pi*theta2/(nx*dx))/2i
b2<-rep(0,2*nx)
c2<-alpha2*nx*dx*base/pi*exp(1i*pi*theta2/(nx*dx))/2i
a1<- -a1
b1<- -b1
c1<- -c1
a2<- -a2
b2<- -b2
c2<- -c2

# #this isn't gonna work...
# #but I think you may be on the right track...
# n<-2*nx
# test.perm<-c(seq(1,n-1,2),seq(2,n,2))
# test.soln<-solve.tridiag.periodic(-dt*a1,1-dt*b1,-dt*c1,posx1[test.perm],2*nx)[test.perm]
# mat<-diag(1-dt*b1)
# diag(mat[-(1:2),])<- -dt*a1[-(1:2)]
# mat[1,n-1]<- -dt*a1[1]
# mat[2,n]<- -dt*a1[2]
# diag(mat[,-(1:2)])<- -dt*c1[-c(n-1,n)]
# mat[n-1,1]<- -dt*c1[n-1]
# mat[n,2]<- -dt*c1[n]
# real.soln<-solve(mat,posx1)
# plot(Re(IFFT(test.soln)))
# plot(Re(IFFT(real.soln)))
# plot(mat%*%test.soln) #it's close...but not quite
# #there's probably something you can do here, but it's not immediately obvious...
# #there actually doesn't really seem to be a way to permute a pentadiagonal matrix into a tridiagonal one...weird!
# #so you'd have to have some kind of "pentadiagonal matrix algorithm" for it to all work...
# #something DEFINITELY for the future time
# #looking into it, I think that the one unstable state at the antioptimum is the way to go
# #just easier and less distortion
# #I can't think of an easy way to introduce the proper "antioptimum" to a quadratic potential...but it might be worth some thought
# #figured it out! Permuting into odd then even order transforms system into two periodic tridiagonal systems
# #which gives a generalizable, divide-and-conquer approach for any period for which n is divisible, if I'm not mistaken...
# #e.g., choosing 4*pi would result in creating 4 periodic tridiagonal matrices by arranging into 4*n, 4*n+1, etc...
# #so might be worth playing around with the more complex trigonometric potential then adding the padded indices to the non-padded ones...
# #but I'm not sure if it's worth it at this point...


foo2<-function(x){
  alpha2*base*(theta2*x+(x[indsll]/12+2*x[inds2]/3-2*x[indsr]/3-x[indsrr]/12)/dbase)
}
#borrowed from a paper...
indsr<-inds1
indsrr<-c(3:(2*nx),1:2)
indsl<-inds2
indsll<-c(2*nx-1,2*nx,1:(2*nx-2))
foo12<-function(x){
  tmpr<-(-3*x+4*x[indsr]-x[indsrr])
  tmpl<-(3*x-4*x[indsl]+x[indsll])
  tmp.inds<-Re(tmpr)<0
  tmpr[tmp.inds]<-tmpl[tmp.inds]
  alpha1*base*(theta1*x-tmpr/(2*dbase))
}
# foo1<-function(x){
#   alpha1*baseish*(theta1*x+(x[inds2]-x[inds1])/(2*dbase))
# }
# foo2<-function(x){
#   alpha2*baseish*(theta2*x+(x[inds2]-x[inds1])/(2*dbase))
# }
BM1<-(-sig21*base2/2)
BM2<-(-sig22*base2/2)

nsteps<-2
t<-1
dt<-t/nsteps
#not worried about efficiency for now
PP<-QQ<-array(Q,c(2,2,2*nx))
QQ[rep(rep(c(TRUE,FALSE),c(1,k)),length.out=k^2)]<-QQ[rep(rep(c(TRUE,FALSE),c(1,k)),length.out=k^2)]+as.vector(rbind(BM1,BM2))
PP[]<-unlist(lapply(asplit(QQ*dt,3),function(ii) t(expm::expm(ii))))
# #not generalized to arbitrary dimensions yet...
# arr.mult<-function(x1,x2){
#   posx1<-PP[1,1,]*x1+PP[1,2,]*x2
#   posx2<-PP[2,1,]*x1+PP[2,2,]*x2
#   list(posx1,posx2)
# }
prex1<-Reduce("+",lapply(1:2,function(ii) BM(1,rnorm(1,0,2),0.1)))
prex1<-prex1/abs(prex1[1])
prex2<-Reduce("+",lapply(1:2,function(ii) BM(1,rnorm(1,0,2),0.1)))
prex2<-prex2/abs(prex2[1])

# #does this help?
# prex1<-Re(IFFT(prex1))
# prex1[prex1<=1e-8]<-1e-8
# prex1<-prex1/sum(prex1)
# prex1<-FFT(prex1)
# prex2<-Re(IFFT(prex2))
# prex2[prex2<=1e-8]<-1e-8
# prex2<-prex2/sum(prex2)
# prex2<-FFT(prex2)
#
# smoother<-BM(1,0,2*sqrt(dt)) #I mean...really helps, but changes dynamics
#definitely not foolproof
for(i in 2:nsteps){
  if(((i-1)%%floor(nsteps/10)==1)){
    plot(Re(fftw::IFFT(T1)[perm.inds])~xx[perm.inds],type="l",ylim=c(0,0.06))
    lines(Re(fftw::IFFT(T3)[perm.inds])~xx[perm.inds],col="red")
    tt<-(i-2)*dt
    tmp<-dnorm(xx,
               15*exp(-alpha*tt)+theta*(1-exp(-alpha*tt)),
               sqrt(sig2/(2*alpha)*(1-exp(-2*alpha*tt))+0.1*exp(-2*alpha*tt)))
    tmp<-tmp/sum(tmp)
    lines(tmp[perm.inds]~xx[perm.inds],col="blue")

    # plot((T2),type="l")
    # lines((BM(1,
    # 15*exp(-alpha*tt)+theta*(1-exp(-alpha*tt)),
    # sig2/(2*alpha)*(1-exp(-2*alpha*tt))+0.1*exp(-2*alpha*tt))),col="blue")
  }

  nsteps<-1000
  t<-1
  dt<-t/nsteps
  T5<-T1<-BM(1,17,0.1)+BM(1,3,0.1)
  const1<- -sig2*base2/2*dt
  const2<-const1-1+alpha*theta*base*dt
  const3<- -alpha*theta*base*dt-1
  i<-1

  #Initial guess
  T2<-T1*BM(dt,0,sig2)
  T2<-T2+foo(T2)*dt

  # T2<-T1*BM(dt,0,sig2)
  # k1<-foo(T2)
  # k2<-foo(T2+dt*k1/3)
  # k3<-foo(T2+dt*(-k1/3+k2/2))
  # k4<-foo(T2+dt*(k1-k2+k3))
  # T2<-T2+dt*(k1+3*k2+3*k3+k4)/8

  #Newton's method
  #the cool thing about this approach is that the number of iterations is related to how inappropriate the timestep is!
  #if it takes a lot of iterations to converge, the timestep is likely too large
  #could use this to actually select a decent timestep for explicit method maybe...
  #might save on computational power overall by eliminating the need to iterate through Newton's method...
  #seems like increasing number of steps by counter/4 to counter/2 works pretty well
  #given a tolerance of 1e-2
  #very nonlinear though...I wonder if there's a more exact way to select a new timstep
  #yeah, counter kinda "blows up" the larger the timestep is
  #it's kind of an exponential decay...but only for sufficient numbers of steps
  #actually, doubling the number of timesteps until you get counter<20 seems like a decent rule...
  nudge<-(T1+T2*const1+foo(T2)*dt-T2)/const2
  counter<-0
  while(max(abs(nudge))>1e-2){
    T2<-T2-nudge/10
    #making sure things don't blow up...
    inds<-abs(T2)>2
    T2[inds]<-2*T2[inds]/abs(T2[inds])
    nudge<-(T1+T2*const1+foo(T2)*dt-T2)/const2
    counter<-counter+1
  }
  # for(j in seq_len(10)){ #just arbitrary for now...should find stopping criterion
  #   #hard because can be unstable and "blow up"
  #   T2<-T2-(T1+T2*const1+foo(T2)*dt-T2)/const2
  #   #this appears to be less stable below...
  #   # T2<-T2-(T1*BM(dt,0,sig2)+foo(T2)*dt-T2)/(alpha*theta*base*dt-1)
  # }
  plot(Re(fftw::IFFT(T2))[perm.inds]~xx[perm.inds],type="l")
  T1<-T2

  # #one last try...
  # #more accurate, but less stable...
  # T3<-T2<-T1*BM(dt,0,sig2)
  # #Newton's method
  # nudge<-(T2+foo(T3)*dt-T3)/(10*const3)
  # while(max(abs(nudge))>1e-4){
  #   T3<-T3-nudge
  #   nudge<-(T2+foo(T3)*dt-T3)/(10*const3)
  # }
  # plot(Re(fftw::IFFT(T3))[perm.inds]~xx[perm.inds],type="l")
  # T1<-T3

  # plot(Re(IFFT(prex1))[perm.inds]~xx[perm.inds],type="l")
  # plot(Re(IFFT(prex1+dt*alpha1*mat%*%prex1))[perm.inds]~xx[perm.inds],type="l")
  # prex1<-prex1+dt*alpha1*mat%*%prex1

  posx1<-PP[1,1,]*prex1+PP[1,2,]*prex2
  # inds11<-abs(posx1)>1e-10
  # posx1[!inds11]<-0
  posx2<-PP[2,1,]*prex1+PP[2,2,]*prex2

  #here we go--fully implicit
  #shit that didn't work...
  #prex1<-mat*prex1+posx1
  #posx1<-prex1-mat*prex1
  #posx1<-(I-mat)*prex1
  #right?
  #oh! just forgot the dt lol
  #slow, but much more stable
  #now the question is how to make this RAPID!
  # prex1<-solve(diag(2*nx)-dt*mat,posx1)
  # prex2<-solve(diag(2*nx)-dt*mat2,posx2)
  #remarkably...it works...really well!
  #even for super long timesteps, it seems to converge to a decent enough solution (at least in the non-padded regions)
  prex1<-solve.tridiag.periodic(-dt*a1,1-dt*b1,-dt*c1,posx1,2*nx)

  # mat<-diag(1-dt*b1)
  # diag(mat[-1,])<- -dt*a1[-1]
  # mat[1,n]<- -dt*a[1]
  # diag(mat[,-1])<- -dt*c1[-n]
  # mat[n,1]<- -dt*c[n]
  # plot(Im(mat%*%prex1)~Im(posx1),type="l")
  # abline(0,1)
  #
  prex2<-solve.tridiag.periodic(-dt*a2,1-dt*b2,-dt*c2,posx2,2*nx)
  #nonperiodic boundaries???
  # prex1<-solve.tridiag(-dt*a1[-1],1-dt*b1,-dt*c1[-2*nx],posx1,2*nx)
  # prex2<-solve.tridiag(-dt*a2[-1],1-dt*b2,-dt*c2[-2*nx],posx2,2*nx)

  # test1<-solve(diag(2*nx)-dt*mat,posx1)
  # test2<-sdetorus::solvePeriodicTridiag(alpha1*base/(2*dbase),1-alpha1*theta1*base,-alpha1*base/(2*dbase),posx1)
  # #damn, doesn't work with complex entries!
  # #so I'd have to do my own implementation...
  # #another time, perhaps...
  # testy<-Matrix::.dense2sparse(mat)
  # #damn, that doesn't work with complex numbers either
  # #hmmm...I have hope for the Thomas algorithm--it doesn't need to be super stable, after all...
  # testy<-solve(diag(2*nx)-dt*mat)

  # tmp<-Re(fftw::IFFT(prex2)[perm.inds])[(nx/2+1):(3*nx/2)]
  # plot(tmp~xx[perm.inds][(nx/2+1):(3*nx/2)],xlim=range(xx[perm.inds][(nx/2+1):(3*nx/2)]),ylim=c(0,max(tmp)),type="l")
  # lines(Re(fftw::IFFT(prex1)[perm.inds])[(nx/2+1):(3*nx/2)]~xx[perm.inds][(nx/2+1):(3*nx/2)],col="red")

  plot(Re(fftw::IFFT(prex2)[perm.inds])~xx[perm.inds],type="l")
  lines(Re(fftw::IFFT(prex1)[perm.inds])~xx[perm.inds],col="red")


  plot(log(Re(IFFT(test))),type="l")
  plot(log(Re(IFFT(posx1))),type="l")

  #holy shit, the solution was so simple potentially?
  #hmmm...helps, but still rather unstable
  #maybe just 1 iteration...
  #still definitely something here, just not as stable as you might hope
  #I think runge kutta's actually making things worse?
  #maybe do a fixed point iteration for each runge kutta step...
#
#   # plot((Re(IFFT(posx1+dt*foo1(posx1)))),type="l")
#   # ppx1<-posx1
#   # ppx1<-posx1+dt*(3*foo1(posx1)/4+foo1(ppx1)/4)
#   # plot((Re(IFFT(ppx1))),type="l")
#   k1<-foo1(posx1)#*inds11
#   k2<-foo1(posx1+dt*k1/3)#*inds11
#   k3<-foo1(posx1+dt*(-k1/3+k2/2))#*inds11
#   k4<-foo1(posx1+dt*(k1-k2+k3))#*inds11
#   prex1<-(posx1+dt*(k1+3*k2+3*k3+k4)/8)#*smoother
#   k42<-foo1(posx1)#*inds11
#   k32<-foo1(posx1-dt*k1/3)#*inds11
#   k22<-foo1(posx1-dt*(-k1/3+k2/2))#*inds11
#   k12<-foo1(posx1-dt*(k1-k2+k3))#*inds11
#   prex1<-(posx1+dt*(k1+3*k2+3*k3+k4)/16+dt*(k12+3*k22+3*k32+k42)/16)#*smoother
#
#   #for x2...
#   k1<-foo2(posx2)#*inds22
#   k2<-foo2(posx2+dt*k1/3)#*inds22
#   k3<-foo2(posx2+dt*(-k1/3+k2/2))#*inds22
#   k4<-foo2(posx2+dt*(k1-k2+k3))#*inds22
#   prex2<-(posx2+dt*(k1+3*k2+3*k3+k4)/8)#*smoother
#   k42<-foo2(prex2)#*inds22
#   k32<-foo2(prex2-dt*k1/3)#*inds22
#   k22<-foo2(prex2-dt*(-k1/3+k2/2))#*inds22
#   k12<-foo2(prex2-dt*(k1-k2+k3))#*inds22
#   prex2<-(posx2+dt*(k1+3*k2+3*k3+k4)/16+dt*(k12+3*k22+3*k32+k42)/16)#*smoother
#
#   plot(Re(fftw::IFFT(prex2)[perm.inds])~xx[perm.inds],type="l")
#   lines(Re(fftw::IFFT(prex1)[perm.inds])~xx[perm.inds],col="red")


  k1<-foo1(posx1)#*inds11
  # plot(Re(IFFT(posx1+dt*foo12(posx1))),type="l")
  # plot(Re(IFFT(posx1+dt*foo1(posx1))),type="l")
  k2<-foo1(posx1+dt*k1/3)#*inds11
  k3<-foo1(posx1+dt*(-k1/3+k2/2))#*inds11
  k4<-foo1(posx1+dt*(k1-k2+k3))#*inds11
  k5<-(k1+3*k2+3*k3+k4)/8

  #you're onto something here...just have to figure out a way to prevent middle low values from "blowing up"
  #looks like a mixture of smoothing and masking indices that appear to have already converged does a decent job...
  #looks like it's mainly the smoothing of derivatives doing the work here...
  #perhaps that's the only way to make this all work
  #so the fixed point iteration stuff is just still unstable, unfortunately...I wonder how I fix that
  # prex1.old<-posx1
  # prex1<-posx1+BM(1,0,0.1)[inds]*dt*k5
  # inds<-rep(TRUE,2*nx)
  # inds<-TRUE
  # inds<-(abs(prex1-prex1.old)>1e-8)&inds
  # prex1.old<-prex1
  # inds<-abs(prex1)>1e-4&inds
  # prex1[inds]<-posx1[inds]+dt*(k5[inds]/3+BM(1,0,0.1)[inds]*2*foo1(prex1)[inds]/3)
  # prex1[inds]<-posx1[inds]+dt*BM(1,0,5*dx)[inds]*foo1(prex1)[inds]#(k5[inds]/2+foo1(prex1)[inds]/2)

  #the below might actually work with some masking of low values...
  prex1<-posx1*BM(dt,0,sig21)
  prex1<-posx1+dt/2*(7*k5/8+foo1(prex1)/8)
  plot(abs(prex1),type="l")
  plot((Re(fftw::IFFT(prex1)[perm.inds]))~xx[perm.inds],type="l")

  k1<-foo2(posx2)
  prex2<-posx2+dt*k1
  prex2<-posx2+dt*(7*k1/8+foo2(prex2)/8)
  plot(abs(prex2),type="l")
  plot(Re(fftw::IFFT(prex2)[perm.inds])~xx[perm.inds],type="l")

  plot(Re(fftw::IFFT(prex2)[perm.inds])~xx[perm.inds],type="l")
  lines(Re(fftw::IFFT(prex1)[perm.inds])~xx[perm.inds],col="red")

  # inds22<-abs(posx2)>1e-10
  # posx2[!inds22]<-0
  #for x1...

  #for implicit solver, must solve equation...
  # finalx1 - foo1(finalx1) = posx1
  # finalx1 - foo1(finalx1) - posx1 = 0
  # finalx1<-posx1+dt*foo1(posx1)
  # posx1+dt*alpha1*mat[-1,-1]%*%prex1=prex1
  # test<-solve(dt*alpha1*mat[-1,-1],prex1[-1]-posx1[-1]) #unfortunately still rather unstable...
  # plot(Re(IFFT(c(1,prex1))),type="l")
  # lines(Re(IFFT(c(1,test))),type="l",col="red") #interesting though...

  #
  # finalx1<-finalx1-2*dbase*(finalx1-dt*foo1(finalx1)-posx1)/(1-dt*alpha1*base*theta1)
  # plot(Re(2*dbase*(finalx1-dt*foo1(finalx1)-posx1)/(1-dt*alpha1*base*theta1)))
  # plot(Re(IFFT(finalx1))) #again you end up with unbounded increase...because you eventually get a negative value...
  # #a matrix-based approach might help...
  # #I feel like there has to be a way to prevent this behavior, but maybe I'm wrong...

  # test<-exp(dt*foo1(posx1))
  # inds<-abs(test)>1
  # testy<-numeric(2*nx)
  # testy[inds]<-posx1[inds]*test[inds]
  # testy[!inds]<-posx1[!inds]+dt*foo1(posx1)[!inds]

  # test<-complex(modulus=abs(posx1)*abs(exp(dt*foo1(posx1))),argument=Arg(posx1)+Arg(dt*foo1(posx1)))
  # plot(Re(IFFT(test)),type="l")
  #
  # lines(Re(IFFT(posx1*exp(dt*foo1(posx1)))),type="l")
  # lines(Re(IFFT(posx1+dt*foo1(posx1))),type="l")
  #
  # plot(abs(1+dt*foo1(posx1)),type="l")
  # plot(abs(exp(dt*foo1(posx1))),type="l")
  # plot(abs(1+dt*foo1(posx1))~abs(exp(dt*foo1(posx1))))
  # abline(0,1)

  k1<-foo1(posx1)#*inds11
  # plot(Re(IFFT(posx1+dt*foo12(posx1))),type="l")
  # plot(Re(IFFT(posx1+dt*foo1(posx1))),type="l")
  k2<-foo1(posx1+dt*k1/3)#*inds11
  k3<-foo1(posx1+dt*(-k1/3+k2/2))#*inds11
  k4<-foo1(posx1+dt*(k1-k2+k3))#*inds11
  prex1<-(posx1+dt*(k1+3*k2+3*k3+k4)/8)#*smoother
  #for x2...
  k1<-foo2(posx2)#*inds22
  k2<-foo2(posx2+dt*k1/3)#*inds22
  k3<-foo2(posx2+dt*(-k1/3+k2/2))#*inds22
  k4<-foo2(posx2+dt*(k1-k2+k3))#*inds22
  prex2<-(posx2+dt*(k1+3*k2+3*k3+k4)/8)#*smoother

  #dividing into two diffusion steps didn't really help...think you really want to use the smoother...
  # prex1<-PP[1,1,]*posx1+PP[1,2,]*posx2
  # prex2<-PP[2,1,]*posx1+PP[2,2,]*posx2

  # plot(abs(prex2),type="l")
  # lines(abs(prex1),col="red")

  # plot(abs(dt*(k1+3*k2+3*k3+k4)/8),type="l")
  # lines(abs(posx2),type="l",col="red")

  # prex1<-Re(IFFT(prex1))
  # prex1[prex1<0]<-0
  # prex1[(nx/2+1):(3*nx/2)]<-0
  # prex1<-FFT(prex1)
  # prex2<-Re(IFFT(prex2))
  # prex2[prex2<0]<-0
  # prex2[(nx/2+1):(3*nx/2)]<-0
  # prex2<-FFT(prex2)

  plot(Re(fftw::IFFT(prex1)[perm.inds])~xx[perm.inds],type="l")
  lines(Re(fftw::IFFT(prex1)[perm.inds])~xx[perm.inds],col="red")

}


#ah, I see the problem! It's not QUITE symmetric because of the rowwise multiplication with base...
#right?
mat[3,]
mat[,3]
#two observations: base can be separated out as a diagonal component...right?
diag(base)%*%mat #etc



nx<-64
dx<-0.1
BM<-BM.DFT(nx,dx)
xx<-c(seq(0,dx*nx,dx),seq(-dx*(nx-1),-dx,dx))
perm.inds<-c((nx+2):(2*nx),1:(nx+1))
base<-.get.base(nx,dx)
base2<-base[[2]]
base<-base[[1]]
dbase<-diff(base[1:2])
inds1<-c(2:(2*nx),1)
inds2<-c(2*nx,1:(2*nx-1))

k<-2
Q<-0.25*matrix(c(-1,1,1,-1),2,2)
alpha1<-0.1
alpha2<-2
theta1<-1
theta2<- -10
sig21<-4
sig22<-1

#matrix form of linear differential operator...
mat2<-mat<-matrix(0,2*nx,2*nx)
diag(mat)<-theta1
mat[cbind(seq_len(2*nx),inds1)]<- -1/(2*dbase)
mat[cbind(seq_len(2*nx),inds2)]<-1/(2*dbase)
diag(mat2)<-base*theta1
mat2[cbind(seq_len(2*nx),inds1)]<- -base/(2*dbase)
mat2[cbind(seq_len(2*nx),inds2)]<- base/(2*dbase)
image(Re(mat))
image(Im(mat))
mat2[2,]
(diag(base)%*%mat)[2,]

tmp<-diag(1/base)-mat
#manual fixing nan's...
tmp[,1]<-c(1,0.5,rep(0,2*nx-3),-0.5)
test<-solve(tmp)
test[,-1]<-test[,-1]%*%diag(1/base[-1])
image(abs(test))
image(abs(solve(diag(2*nx)-mat2)))
#ugh, it almost works, but the 1/0 mucks everything up!
#also, still don't really know what to do about solving mat in the first place...
#but there's definitely some patterns that may be used it feels like

eig<-eigen(mat)
plot(eig$values,type="l") #ooh, all real, sine curve?
image(Im(eig$vectors)) #look like all sine waves
#so the issue is that this gets...wonky if multiplied by diagonal matrix...
plot((eig$vectors[,12]),type="l") #not quite sine waves...they're weird!
#odd indices are same as even indices just offset
#seem to be sine waves of certain frequency that decay periodically
#phase seems kinda random, but periodicity always seems to be 4 here...weird!
image(abs(solve(eig$vectors)))
image(abs(t(eig$vectors))) #does seem to be inverse to its conjugate transpose
image(Re(eig$vectors%*%t(Conj(eig$vectors)))) #sweet

#hmmm...I think there is a connection here but difficult to see
#hard because eigenvalues don't correspond to particular rows...
plot(sort(Re(eig$values)))
eig2<-eigen(mat2)
plot(sort(Im(eig2$values)))
image(abs(eig2$vectors))
image(abs(eig$vectors))

#well this identity at least works...
image(abs(mat[-1,-1]))
image(abs(solve(mat[-1,-1])))
image(abs(mat2[-1,-1]))
image(abs(solve(mat2[-1,-1])))
image(abs(solve(mat[-1,-1])%*%diag(1/base[-1])))
#so theoretically, this should be the inverse...
tmp<-solve(mat)%*%diag(1/base)
image(abs(solve(mat)%*%diag(1/base)))
plot(Re(tmp[2,]),type="l") #probs meant to be 1s in first column?
tmp[,1]<-1
image(abs(mat2%*%tmp)) #nah...
tmp[,1]<-0
image(abs(mat2%*%tmp)) #nearly so...
#probably technically meant to be 1/0 in the first entry of the first column...
#not sure how this helps with the matrix exponential per se...but interesting!
#unfortunately doesn't seem to help too too much with eigendecomp--you get some interesting inequalities but no identities...


image(Re((eig$vectors%*%diag(exp(0.1*eig$values))%*%t(Conj(eig$vectors)))*base)) #doesn't really seem to do much...
testy<-Re(eig2$vectors%*%diag(exp(0.1*eig2$values))%*%(solve(eig2$vectors)))
image(Re(eig2$vectors%*%diag(exp(0.1*eig2$values))%*%(solve(eig2$vectors))))
plot(Re(fftw::IFFT(testy[,64])),type="l")
plot(sin(seq(0,99*2*pi,length.out=129)[-129]),type="l")
image(abs(eig$vectors))
image(abs(eig2$vectors))

image(Re(eig2$vectors))
image(Re(t(eig2$vectors)))
testy<-solve(eig2$vectors)
image(Re(solve(eig2$vectors)))
image(Re(t(eig2$vectors)%*%diag(1/base)))
Re(t(eig2$vectors)%*%diag(1/base))
#feels like a deadend, but an interesting investigation nonetheless...
#there's definitely some patterns to it all, but it seems to complicated to really make sense of

#wait second...
tmp<-solve(mat,eig$vectors)
image(abs(tmp)) #interesting...
#so you can find a matrix that transforms the base matrix into its eigenvectors
#which should immediately give a recipe for forming the appropriate eigendecompostion
#the question is if "tmp" above has a predictable form...
plot(Re(tmp[,100]),type="l") #huh, maybe? certainly not intuitive, but might be something there
image(abs(diag(base)%*%mat%*%tmp)) #whoa...



#imaginary diagonal
#real off-diagonal, negatives of each other
test<-eigen(mat)
image(Im(test$vectors))
plot(rev(Im(test$vectors[,5])))
plot(-(Im(test$vectors[,6])))
#well, first clue is that each even eigenvector is a reflection of the complex conjugate of the last
plot(Im(test$values))
plot(diag(mat))
#I think it's simply two straight lines, starting with negative...
#oh not quite--because setting theta to 0 destroys everything...
#might actually be more like some kind of exponential decay from 2*base[2*nx+1]/dbase
#though this would make it constant across values of theta...which it decidedly isn't
#it actually might be square root curves...
lines(-64*sqrt((1:128)/128)+64)
# test.eigvals<-rep(diag(mat)[(nx+1):(2*nx)],each=2)*c(-1,1)
#Now the question is what the HECK the eigenfunctions are! Definitely a pattern, but not clear
plot(abs(test$vectors[,100]),type="l") #definitely a consistent shape...
plot(abs(fftw::IFFT(test$vectors[,50]))[perm.inds]~xx[perm.inds],type="l") #such a weird, hard-to-characterize shape though...
abline(v=0)

testy<-solve(test$vectors)
image(Re(testy))
image(Re(t(Conj(test$vectors))))
image(Re(test$vectors%*%testy))
image(Re(test$vectors%*%t(Conj(test$vectors)))) #weird...that's not a numerical error is it?

#good news is that the eigenvectors are indeed the transpose of their own inverse...
#oh they should
#which saves on computation for sure
#but it still doesn't answer how to compute the eigenfunctions...
#damn, not true in all cases again, as you were starting to suspect...
test.result<-test$vectors%*%diag(exp(1*test$values))%*%solve(test$vectors)
image(abs(test.result))
test.result2<-test$vectors%*%diag(exp(0.1*test.eigvals))%*%t(test$vectors)
image(abs(test.result2)) #not quite...but could be numerical errors to do with the eigenvectors
#so what are the eigenvectors *trying* to be?
plot((abs(test$vectors[,8])),type="l")
plot(Im(fftw::IFFT(test$vectors[,5]))[perm.inds]~xx[perm.inds],type="l")
#it's interesting...there's definitely a pattern, but I'm unsure what it is precisely...

plot(abs(test.result[,7]),type="l")
plot(Re(fftw::IFFT(test.result[,50]))[perm.inds]~xx[perm.inds],type="l") #seems to be sinusoidal wave with amplitude of 1/(2*nx) and frequency of column-1
abline(v=2.5)
#Im is the sine part
#Re is the cosine part
#pattern becomes way less clear as dt gets large though...
#seems like it starts to reach some kind of stationary distribution....
#kinda become sinc-like...
#interesting...
image(abs(test$vectors))
plot(abs(test$vectors[,3]),type="l")

