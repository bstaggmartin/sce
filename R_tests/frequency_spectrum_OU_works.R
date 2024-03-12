#1/10: alright, so here's the plan:
#use an implicit "test step" and use the number of iterations to determine whether the current step size must be increased or decreased
#if the number of iterations is too high (>30 seems like a decent rule of thumb for a tolerance of 1e-2), divide step size in half
#if the number of iterations is too low (<20 maybe?), double the step size
#there's probably a less crazy way to do this, but this seems like a nice, dynamic way to adjust the step size on the fly and use the explicit scheme
#overall reduces the number of evaluations and prevents undesirable behavior...
#I think I like this method because it's more straight-forward, more accurate, and let's the particular initial conditions speak for themselves, rather than trying to do anything super fancy with the math...
#Yeah, super general purpose too...as long as you can find good cutoffs for the counters
#Oh, you're also dividing the nudges by 10 for additional stability!
#Damn, quite a few layers deep...but I think this is the way to go
#Will generally save time, or at least memory in case the number of timesteps exceeds 2*nx...
#Definitely saves on time if you have state-dependent OU...
#Also allows you to use matrix exponential trick since it's explicit!


####RANDOM OU TEST####
#holy crap it actually works!
nx<-1024
dx<-0.1
BM<-BM.DFT(nx,dx)
xx<-c(seq(0,dx*nx,dx),seq(-dx*(nx-1),-dx,dx))
perm.inds<-c((nx+2):(2*nx),1:(nx+1))


#taking derivatives...
base<-.get.base(1024,0.1)[[1]]
plot(Re(fftw::IFFT(BM(10.025,0,1)))[perm.inds]~xx[perm.inds],type="l")
lines(Re(fftw::IFFT(BM(10.025,0,1)*(-base)))[perm.inds]~xx[perm.inds],type="l")
lines(Re(fftw::IFFT(BM(10.025,0,1)*(-base)^2))[perm.inds]~xx[perm.inds])
lines(Re(fftw::IFFT(BM(10.025,0,1)*(-base)^3))[perm.inds]~xx[perm.inds])
#suggests a simple pseudospectral recipe, right?
#not as much as I thought...but I found a recipe that works
#not sure about numerical stability, though, and doesn't preserve integral of function
#figured out potentially a way to avoid repeatedly fft'ing, but still working on it
#seems highly numerically unstable and can't figure out how to shift optimum to non-zero value...
# pot<-dnorm(xx,15,3)[-((nx/2+1):(3*nx/2))]
# plot(pot[perm.inds]~xx[perm.inds])
# pot<-pot/sum(pot)
nsteps<-1200
t<-30
dt<-t/nsteps
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
#A potential hint: bubbles always seem to occur around x=20 with dx=0.05 and nx=1024
#hmmm...this is essentially numerical integration--maybe adopting some runge-kutta techniques would help deal with the numerical instability...
for(i in 2:nsteps){
  T2<-BM(dt,0,2)*T1
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

#oof, explicit runge kutta is pretty bad...
#wow, this actually works really pretty well now!
#so you want to take difference in amount of diffusion into account, but how exactly is an open question...
#the main takeaway is that some kind of runge kutta technique is a good idea
#seems pretty stable now!
#can do some precomputation too for extra speed...done
#so now the question is how to determine the optimal step size to prevent the artifacts?
#yeah, much better, but how to know when it's going to blow up? Must be some way to find optimal step size...
#alright:
# - only doing runge kutta on the non-diffusion step seems more stable
# - but still blows up under mysterious conditions...

theta<- -10
alpha<-1
sig2<-10

nsteps<-201
t<-10
dt<-t/nsteps
dbase<-diff(base[1:2])
inds1<-c(2:(2*nx),1)
inds2<-c(2*nx,1:(2*nx-1))
foo<-function(x){
  -alpha*base*((x[inds1]-x[inds2])/(2*dbase)-theta*x)
}
full.kern<-BM(dt,0,sig2)

T1<-BM(1,15,0.1)
#requires similar time steps to that of classic method
for(i in 2:nsteps){
  if(((i-1)%%floor(nsteps/10)==1)){
    plot(Re(fftw::IFFT(T1)[perm.inds])~xx[perm.inds],type="l")
  }
  T2<-T1*full.kern
  k1<-foo(T2)
  k2<-foo(T2+0.4*dt*k1)
  k3<-foo(T2+dt*(0.29697761*k1+0.15875964*k2))
  k4<-foo(T2+dt*(0.21810040*k1-3.05096516*k2+3.83286476*k3))
  T1<-T2+dt*(0.17476028*k1-0.55148066*k2+1.20553560*k3+0.17118478*k4)
  # plot(Re(fftw::IFFT(T1)[perm.inds])~xx[perm.inds],type="l")
}

T1<-BM(1,15,0.1)
half.kern<-BM(dt/2,0,sig2)
for(i in 2:nsteps){
  if(((i-1)%%floor(nsteps/10)==1)){
    plot(Re(fftw::IFFT(T1)[perm.inds])~xx[perm.inds],type="l")
  }
  # T2<-T1*half.kern
  T3<-T1*full.kern
  k1<-foo(T3) #in a true runge kutta method, this would be T1...but that actually makes it less stable...
  #so should you be conservative with diffusion part?
  #that actually does greatly improve things...interesting!
  k2<-foo(T3+dt*k1/2)
  k3<-foo(T3+dt*k2/2)
  k4<-foo(T3+dt*k3)
  T1<-T3+dt*(k1+2*k2+2*k3+k4)/6
  # plot(Re(fftw::IFFT(T1)[perm.inds])~xx[perm.inds],type="l")
}

T2<-T1<-BM(1,15,0.1)
foo2<-function(dt,x){
  tmp1<-x*BM(dt,0,sig2)
  dt*foo(tmp1)+tmp1-x
}
base2<- -base^2
for(i in 2:nsteps){
  if(((i-1)%%floor(nsteps/10)==1)){
    plot(Re(fftw::IFFT(T1)[perm.inds])~xx[perm.inds],type="l")
    lines(Re(fftw::IFFT(T2)[perm.inds])~xx[perm.inds],col="red")
    tt<-(i-2)*dt
    tmp<-dnorm(xx,
               15*exp(-alpha*tt)+theta*(1-exp(-alpha*tt)),
               sqrt(sig2/(2*alpha)*(1-exp(-2*alpha*tt))+2*exp(-2*alpha*tt)))
    tmp<-tmp/sum(tmp)
    lines(tmp[perm.inds]~xx[perm.inds],col="blue")

    # plot((T2),type="l")
    # lines((BM(1,
            # 15*exp(-alpha*tt)+theta*(1-exp(-alpha*tt)),
            # sig2/(2*alpha)*(1-exp(-2*alpha*tt))+0.1*exp(-2*alpha*tt))),col="blue")
  }
  #this is perhaps more stable, but not as accurate...
  k1<-foo2(dt,T1)
  k2<-foo2(dt/2,T1+k1/2)
  k3<-foo2(dt/2,T1+k2/2)
  k4<-foo2(dt,T1+k3)
  T1<-T1+(k1+2*k2+2*k3+k4)/6

  # let's try more explicitly differentiating the convolution step...
  plot(abs(T2))
  plot(abs(T2*exp(-sig2*base2/2*dt)))
  plot(Im((T2*exp(-sig2*base2/2*dt)-T2)/dt))
  plot(abs(fftw::IFFT(T2*full.kern)),type="l")
  plot(abs(fftw::IFFT(T2+T2*-sig2*base2/2*dt)),type="l")

  #Newton's method for diffusion...
  T3<-T2
  T3<-T3-(T2+T3*-sig2*base2/2*dt-T3)/(-sig2*base2/2*dt-1) #wow, Newton's method works really well here...
  plot(abs(T3))

  #now let's try adding in the OU component...
  #this might actually work alright...
  #seems more stable, but not perfect--one more try
  T3<-T2#+T2*sig2*base2/2*dt+foo(T2)*dt
  T3<-T3-(T2+T3*-sig2*base2/2*dt+foo(T3)*dt-T3)/(-sig2*base2/2*dt-1+2*alpha*theta*base*dt)
  #this made it less stable...
  # T3<-T3-
  #   (T2+T3*-sig2*base2/2*dt+foo(T2+T3*-sig2*base2/2*dt)*dt-T3)/
  #   (-sig2*base2/2*dt-1+alpha*theta*sig2*base^3*dt^2)
  plot(abs(T3))
  plot(Re(fftw::IFFT(T3))[perm.inds]~xx[perm.inds],type="l")

  T3<-T2+T3*-sig2*base2/2*dt
  plot(abs(T3),type="l")
  plot(abs(fftw::IFFT(T3)),type="l")
  plot(Im(T2*-sig2*base2/2))




  T3<-T2
  T3<-T2*full.kern
  k1<-foo(T3)
  k2<-foo(T3+dt*k1/3)
  k3<-foo(T3+dt*(-k1/3+k2/2))
  k4<-foo(T3+dt*(k1-k2+k3))
  T2<-T3+dt*(k1+3*k2+3*k3+k4)/8


  # plot(Re(fftw::IFFT(T1)[perm.inds])~xx[perm.inds],type="l")
}

#this is working pretty well, but the alpha's pull feels just a tad too strong...
#fixed! It now seems to actually work...
#will be trickier to extend to state-dependent domain but seems like a good starting point
#still have no idea how to judge when stable or not, but seems to allow for much larger timesteps...
####THIS IS THE BEST VERSION, USING IMPLICIT SCHEME WITH TRUNCATED NEWTON'S METHOD####
T5<-T1<-BM(1,15,5)+BM(1,-5,0.2)
const1<- -sig2*base2/2*dt
# #higher order approximations don't seem to do much...
# const1<-(const1+const1^2/2+const1^3/6+const1^4/24+const1^5/120)*dt
#intermediate wide normal seems to result from time discretization...
#fascinating--I actually think it's from the circular boundary conditions potentially!
#no, circular boundary conditions are really only a problem when x0 or theta is beyond the central half, as far as I can discern...
#circular boundary conditions definitely start to have a pretty drastic effect with the backwards process...
#on the bright side, the backwards process is a simple as flipping the sign of foo
#you'll just have to think more on how to prevent the weird circular influences...
#additional zero-padding seems to do the trick and also improves stability overall...
#which begs the question of how to find the best counter cutoff...
#oh no, just a typo! counter cutoff around 30 still works...
#though it undoubtedly still depends on a myriad of factors
#another thing to recognize is that the oscillations tend to be associated with the absolute value exceeding 1...
#doesn't prevent oscillations, but does seem to keep them bounded...
const2<-const1-1+alpha*theta*base*dt
const3<-alpha*theta*base*dt-1
i<-1
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

  dbase<-diff(base[1:2])
  inds1<-c(2:(2*nx),1)
  inds2<-c(2*nx,1:(2*nx-1))
  foo<-function(x){
    -alpha*base*((x[inds1]-x[inds2])/(2*dbase)-theta*x)
  }
  base2<- -base^2

  nsteps<-100
  t<-10
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


  T4<-T5*BM(dt,0,sig2)
  k1<-foo(T4)
  k2<-foo(T4+dt*k1/3)
  k3<-foo(T4+dt*(-k1/3+k2/2))
  k4<-foo(T4+dt*(k1-k2+k3))
  T5<-T4+dt*(k1+3*k2+3*k3+k4)/8

  plot(Re(fftw::IFFT(T1)[perm.inds])~xx[perm.inds],type="n",ylim=c(0,0.06))
  lines(Re(fftw::IFFT(T5)[perm.inds])~xx[perm.inds],col="red")
  abline(v=10)
  i<-i+1
  tt<-(i-1)*dt
  # tmp<-dnorm(xx,
  #            3*exp(-alpha*tt)+theta*(1-exp(-alpha*tt)),
  #            sqrt(sig2/(2*alpha)*(1-exp(-2*alpha*tt))+0.1*exp(-2*alpha*tt)))+
  #   dnorm(xx,
  #         17*exp(-alpha*tt)+theta*(1-exp(-alpha*tt)),
  #         sqrt(sig2/(2*alpha)*(1-exp(-2*alpha*tt))+0.1*exp(-2*alpha*tt)))
  tmp<-dnorm(xx,
             3*exp(alpha*tt)+theta*(1-exp(alpha*tt)),
             sqrt(sig2*(exp(2*alpha*tt)-1)/(2*alpha)+0.1*exp(2*alpha*tt)))+
    dnorm(xx,
          17*exp(alpha*tt)+theta*(1-exp(alpha*tt)),
          sqrt(sig2*(exp(2*alpha*tt)-1)/(2*alpha)+0.1*exp(2*alpha*tt)))
  tmp<-2*tmp/sum(tmp)
  lines(tmp[perm.inds]~xx[perm.inds],col="blue")
}

out<-7:500
for(k in seq_along(out)){
  nsteps<-10*(k+6)
  t<-10
  dt<-t/nsteps
  T5<-T1<-BM(1,15,0.1)+BM(1,-5,0.1)
  const1<- -sig2*base2/2*dt
  const2<-const1-1+alpha*theta*base*dt
  const3<-alpha*theta*base*dt-1
  i<-1

  #Initial guess
  T2<-T1*BM(dt,0,sig2)
  T2<-T2+foo(T2)*dt
  nudge<-(T1+T2*const1+foo(T2)*dt-T2)/const2
  counter<-0
  while(max(abs(nudge))>1e-2){
    T2<-T2-nudge/10
    nudge<-(T1+T2*const1+foo(T2)*dt-T2)/const2
    counter<-counter+1
  }
  out[k]<-counter
}


nsteps<-1001
t<-10
dt<-t/nsteps
#slightly better...
#can get away with fewer timesteps...
T1<-BM(1,15,0.1)
full.kern<-BM(dt,0,sig2)
for(i in 2:nsteps){
  if(((i-1)%%floor(nsteps/10)==1)){
    plot(Re(fftw::IFFT(T1)[perm.inds])~xx[perm.inds],type="l",ylim=c(0,0.06))
  }
  T2<-T1*full.kern
  k1<-foo(T2)
  k2<-foo(T2+dt*k1/3)
  k3<-foo(T2+dt*(-k1/3+k2/2))
  k4<-foo(T2+dt*(k1-k2+k3))
  T1<-T2+dt*(k1+3*k2+3*k3+k4)/8
  # plot(Re(fftw::IFFT(T1)[perm.inds])~xx[perm.inds],type="l")
}

#backwards euler maybe?
#not stable--pde is too stiff
#probably something here, but it will take some tinkering and thinking to derive gradients and such...
#also does optim work with complex numbers? probs not...
#goes haywire eventually, but I think I'm onto something with the Newtony method here
fooy<-function(par){
  sum(abs(T1*full.kern+dt*foo(par)-par))
}
T1<-BM(1,15,0.1)
for(i in 2:nsteps){
  if(((i-1)%%floor(nsteps/10)==1)){
    plot(Re(fftw::IFFT(T1)[perm.inds])~xx[perm.inds],type="l")
  }

  #perhaps matrix-based?
  T2<-T1*full.kern
  mat<-matrix(0,2*nx,2*nx)
  diag(mat)<-1*-2*base*10
  mat[cbind(1:(2*nx),inds2)]<-1*base*10/(2*dbase)
  mat[cbind(1:(2*nx),inds1)]<- -1*base*10/(2*dbase)
  test<-solve(mat[-1,-1],(T1*full.kern+dt*foo(T2)-T2)[-1])
  ran<-range(which(abs(Re(test))<1e-4))
  test[ran[1]:ran[2]]<-0

  plot(Re(test))
  T2[-1]<-T2[-1]-test #this actually works rather well, but ends up running into OTHER artifacts...
  plot(Re(fftw::IFFT(T2)[perm.inds])~xx[perm.inds],type="l")
  lines(Re(fftw::IFFT(T1*full.kern+dt*foo(T1*full.kern)))[perm.inds]~xx[perm.inds],type="l",col="red")
  T2<-T1

  T2.test<-T2
  #pshhh, I can't figure this one out



  #I think the problem is this simply doesn't make sense...
  #You will never find a perfect fit, and numerical errors accumulate...
  #It does do...something, but not sure it's the correct something
  T2<-T1*full.kern
  tmp<-(T1*full.kern+dt*foo(T2)-T2)/(1*-2*10*base)
  tmp[is.nan(tmp)]<-0 #this is the problem--it's not technically zero for some reason...
  #the thing I don't understand is that this really SHOULD stay 1, right?
  #so things integrate to 1?
  #the key is to probably nudge things such that the derivative is roughly linear in that area...
  # tmp[is.nan(tmp)]<-abs(tmp[2])-(abs(tmp[3])-abs(tmp[2]))
  #hmmm...I think the trick is to rescale to make sure things don't blow up
  # plot(abs(tmp))
  T2<-T2-tmp #don't know why it's plus instead of minus...
  T2[-1]<-T2[-1]/T2[2]
  # T2<-T2/T2[1]
  # plot(abs(T2))
  # tmp<-(T1*full.kern+dt*foo(T2)-T2)/(1*-2*10*base)
  # T2<-T2-tmp
  plot(Re(fftw::IFFT(T2)[perm.inds])~xx[perm.inds],type="l")
  lines(Re(fftw::IFFT(T1*full.kern+dt*foo(T1*full.kern)))[perm.inds]~xx[perm.inds],type="l",col="red")

  #It actually does seem to work, but you need to figure out a good stopping criterion--it eventually runs into problems
  #Probably something with the rescaling step that's more stable...

  # tmp<-T1*full.kern+dt*foo(T2)-T2
  # dtmp<-(foo(T2+1e-4i)-foo(T2))/1e-4i
  # plot(-dt*2*10*Im(base)) #I believe the derivative operation kinda cancels out...if you treat everything on a case-by-case basis at least
  # plot(dt*Im(dtmp))

  test<-optim(T1,fooy)
  T2<-T1*full.kern+dt*foo(T2)
  plot(Re(fftw::IFFT(T2)[perm.inds])~xx[perm.inds],type="l")
  T1<-T2


  k1<-foo(T2)
  k2<-foo(T2+0.4*dt*k1)
  k3<-foo(T2+dt*(0.29697761*k1+0.15875964*k2))
  k4<-foo(T2+dt*(0.21810040*k1-3.05096516*k2+3.83286476*k3))
  T1<-T2+dt*(0.17476028*k1-0.55148066*k2+1.20553560*k3+0.17118478*k4)
  # plot(Re(fftw::IFFT(T1)[perm.inds])~xx[perm.inds],type="l")
}

#maybe isolate the runge kutta step to the non-diffusion part?
#this works slightly better, but something still seems off...
#seems to work far from optimum value, but not close to it...
#figured it out! typo
#reduces number of required steps by ~1/3
nsteps<-900
t<-10
dt<-t/nsteps
T1<-BM(1,15,0.1)
dbase<-diff(base[1:2])
inds1<-c(2:(2*nx),1)
inds2<-c(2*nx,1:(2*nx-1))
for(i in 2:nsteps){
  T2<-BM(dt,0,2)*T1
  dT2<-(T2[inds1]-T2[inds2])/(4*dbase)
  k1<-1*((dT2+10*T2)*-2*base)
  tmp<-T2+k1*dt/2
  dtmp<-(tmp[inds1]-tmp[inds2])/(4*dbase)
  k2<-1*((dtmp+10*tmp)*-2*base)
  tmp<-T2+k2*dt/2
  dtmp<-(tmp[inds1]-tmp[inds2])/(4*dbase)
  k3<-1*((dtmp+10*tmp)*-2*base)
  tmp<-T2+k3*dt
  dtmp<-(tmp[inds1]-tmp[inds2])/(4*dbase)
  k4<-1*((dtmp+10*tmp)*-2*base)
  T1<-T2+dt/6*(k1+2*k2+2*k3+k4)
  # T1<-T2+dt*k1
  # plot(Re(fftw::IFFT(T1)[perm.inds])~xx[perm.inds],type="l")
  if(i%%100==0){
    plot(Re(fftw::IFFT(T1)[perm.inds])~xx[perm.inds],type="l")
  }
}


nsteps<-101
t<-10
dt<-t/nsteps


T5<-T1<-BM(1,15,5)+BM(1,-5,0.2)
const1<- -sig2*base2/2*dt
# #higher order approximations don't seem to do much...
# const1<-(const1+const1^2/2+const1^3/6+const1^4/24+const1^5/120)*dt
#intermediate wide normal seems to result from time discretization...
#fascinating--I actually think it's from the circular boundary conditions potentially!
#no, circular boundary conditions are really only a problem when x0 or theta is beyond the central half, as far as I can discern...
const2<-const1-1+alpha*theta*base*dt
const3<-alpha*theta*base*dt-1
i<-1
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


  #Initial guess
  T2<-T1*BM(dt,0,sig2)
  T2<-T2+foo(T2)*dt
  #Newton's method
  for(j in seq_len(10)){ #just arbitrary for now...should find stopping criterion
    #hard because can be unstable and "blow up"
    T2<-T2-(T1+T2*const1+foo(T2)*dt-T2)/const2
  }
  T1<-T2

  plot(Re(fftw::IFFT(T1)[perm.inds])~xx[perm.inds],type="l")
  lines(Re(fftw::IFFT(T4)[perm.inds])~xx[perm.inds],col="red")
}
