#6/20: after a VERY frustrating day, I have finally figured out the problem
#annoyingly enough, you didn't have to update R or anything
#there seems to be some kind of memory allocation error associated with fftw's DCT code if plan isn't appropriately specified ahead of time...
#I can't see why this would be happening though...

#6/24: thought to mix dcts with dsts to get desired results...but it doesn't work
#dst convolutions are extremely odd to wrap your head around--they exhibit odd symmetry like intended, but they still yield symmetric convolutions
#instead of getting the pos/neg copies of each conv. kernel moving in different directions, you get two mirror copies of pos/neg kernel moving in both directions
#(which is neat, but not what you're looking for)

#it was good to experiment with this and wrap your head around it all, but I think padded DFTs are unfortunately your best bet--there really isn't another way
#DCT's encode reflective boundaries way more naturally, but they're just too limited in their ability to do asymmetric convolutions
#And this is a big issue for your OU stuff, because you need the asymmetry to numerically calculate derivatives and such
#Could, if you drop the drift components, use for SCE, but feels unnecessary

install.packages("gsignal")
library(gsignal)

nx<-1024
dx<-0.01
xx<-seq(0,by=dx,length.out=nx)
base<-seq(0,pi/dx,length.out=nx)
BM.DCT<-function(t,rate){
  exp(-t*rate/2*base^2)
}
seed<-dnorm(xx,5,0.1)
seed[seed<1e-15]<-1e-15

dct.mat<-rbind(0.5,
               cos(pi/(nx-1)*outer(seq_len(nx-2),0:(nx-1))),
               c(0.5,-0.5))
DCT<-function(x){
  dct.mat%*%x
}
IDCT<-function(x){
  2/(nx-1)*dct.mat%*%x
}
plot(DCT(seed))

plot(-DCT(c(0,DST(seed),0)),type="l")
plot(diff(seed),type="l")

plot(IDCT(DCT(seed)*BM.DCT(1,1))~xx) #sweet
dst.mat<-sin(pi/(nx+1)*outer(seq_len(nx-2),seq_len(nx-2)))
DST<-function(x){
  dst.mat%*%x[-c(1,nx)]
}
IDST<-function(x){
  c(0,2/(nx+1)*dst.mat%*%x,0)
}
plot(DST(seed))
plot(IDST(DST(seed)))
plot(IDST(DST(seed)*BM.DCT(1,1)[-c(1,nx-1)])~xx)

plot(IDCT(DCT(seed)*BM.DCT(1,0.2)*cos(2*base))~xx)
t<-1
plot(IDST(DST(seed)*BM.DCT(t,0.2)[-c(1,nx-1)]*sin(t*2*base[-c(1,nx-1)]))~xx)
plot(IDST(DST(seed)*BM.DCT(t,0.2)[-c(1,nx-1)]*sin(t*2*base[-c(1,nx-1)]))+
       IDCT(DCT(seed)*BM.DCT(1,0.2)*cos(2*base))~xx)

plot(IDCT(DCT(seed)*DCT(dnorm(xx,2,0.1)))~xx) #boundaries messed up...but shape is there
plot(IDST(DST(seed)*DST(dnorm(xx,0.3,0.3)))~xx)
plot(IDST(DST(dnorm(xx,2,0.1))))

dst.mat<-rbind(sin(pi/nx*outer(seq_len(nx-3),seq(0.5,nx-2.5))),
               c(0.5,-0.5))
plot(DST(seed))
idst.mat<-sin(pi/nx*outer(seq(0.5,nx-1.5),seq_len(nx-2)))
IDST<-function(x){
  c(0,2/nx*idst.mat%*%x,0)
}
plot(IDST(DST(seed)*DST(dnorm(xx,3,0.3))))

plan<-planDCT(n=nx,type=1)
test<-DCT(seed,type=1,plan=plan)

Im(exp(2*1i*base))[2]
sin(2*base)[2]

library(fftw)

nx<-1024
dx<-0.01
xx<-seq(0,by=dx,length.out=nx)
base<-seq(0,pi/dx,length.out=nx)
BM.DCT<-function(t,rate){
  exp(-t*rate/2*base^2)
}

seed<-dnorm(xx,5,0.1)
seed[seed<1e-15]<-1e-15
#the DCT function appears to be bugged for some reason!
plan<-planDCT(n=nx,type=1)
test<-DCT(seed,type=1,plan=plan)
plot(log(IDCT(test*BM.DCT(1,1),plan=plan))~xx,type="l") #sweeet, reflective boundaries!

#onto something here...but doesn't quite cancel out
#must be missing something somewhere
#but this may be on the right track somehow...
tmp<-test*BM.DCT(1,1)
tmp<-tmp[c(1024,1:1023)]*Im(exp(3*1i*base))
plot(IDCT(test*BM.DCT(1,1)*Re(exp(3*1i*base))-tmp,plan=plan)/2~xx)

tmp<-c(0,(test*BM.DCT(1,1))[-1024])*Im(exp(2*1i*base))
plot(IDCT(test*BM.DCT(1,1)*Re(exp(2*1i*base))-tmp,plan=plan)~xx)

#yeah, can't figure out if asym kernels can work...
#feels like there could be something here, but not that I can figure out
base2<-1i*seq(0,-2*pi/dx,length.out=2*nx+1)[-(2*nx+1)]

t<-0

plot(IDCT(test*BM.DCT(t,1)*Re(exp(t*base2)),plan=plan)~xx)
t<-t+1



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
