#I think all scaling is correct now, but I might've messed up somewhere

.get.base<-function(nx,dx){
  base<-seq(0,pi/dx,length.out=nx+1)^2
  c(base,rev(base)[-c(1,nx+1)])
}

#Brownian Motion process DFT
#t is time, sig2 is variance parameter
#' @export
BM.DFT<-function(nx,dx){
  base<- -0.5*.get.base(nx,dx)
  out<-function(t,sig2){
    exp(t*sig2*base)
  }
}

#compound possion process DFT
#t is time, rate is lambda, jump is variance of normal distribution
#' @export
JN.DFT<-function(nx,dx){
  base<- -0.5*.get.base(nx,dx)
  out<-function(t,rate,jump){
    exp(t*rate*(exp(jump*base)-1))
  }
  out
}

#something broken here
#ahhh, typo--instead of alpha^2, it should be alpha on pg. 3 of landis and schraiber supplement
#now seems right
#normal inverse gaussian process DFT
#t is time, rate is delta, jump is 1/alpha2
#' @export
NIG.DFT<-function(nx,dx){
  base<-.get.base(nx,dx)
  out<-function(t,rate,jump){
    inv.jump<-1/jump
    exp(t*rate*(sqrt(inv.jump)-sqrt(inv.jump+base)))
  }
  out
}

#this seems right
#variance gamma process DFT
#t is time, rate is 1/nu, jump is tau2
#' @export
VG.DFT<-function(nx,dx){
  base<-.get.base(nx,dx)
  out<-function(t,rate,jump){
    (1+base*jump/(2*rate))^(-t*rate)
    #exp(t*x) = (1+base*jump/(2*rate))^(-t*rate)
    #t*x = -t*rate*ln(1+base*jump/(2*rate))
    #x = -rate*ln(1+base*jump/(2*rate))
  }
  out
}

#dirac delta DFT
#x is the location parameter, x0 is the minimum value of the considered range
#' @export
dirac.DFT<-function(nx,dx,x0=0){
  len<-2*nx+1
  base<- -1i*seq(0,2*pi/dx,length.out=len)[-len]
  out<-function(x){
    exp((x-x0)*base)
  }
}

#normal distribution DFT
#x is location parameter, sig2 is variance parameter, x0 is the minimum value of the considered range
#equivalent to convolution of dirac delta and BM DFTs
#' @export
gauss.DFT<-function(nx,dx,x0=0){
  tmp<-seq(0,pi/dx,length.out=nx+1)
  loc.base<- -1i*c(tmp,tmp[-c(1,nx+1)]+tmp[nx+1])
  tmp<-tmp^2
  scale.base<- -0.5*c(tmp,rev(tmp)[-c(1,nx+1)])
  out<-function(x,sig2){
    exp((x-x0)*loc.base+sig2*scale.base)
  }
}
#
# nx<-1024
# dx<-0.245
# xx<-seq(-nx*dx/2,nx*dx/2-dx,length.out=nx)
#
# tmp<-seq(0,pi/dx,length.out=nx/2+1)
# loc.base<- -1i*c(tmp,tmp[-c(1,nx/2+1)]+tmp[nx/2+1])
# tmp<-tmp^2
# scale.base<- -0.5*c(tmp,rev(tmp)[-c(1,nx/2+1)])
# plot(fftw::IFFT(exp((0-xx[1])*loc.base+1*scale.base))~xx)

# nx<-1024
# dx<-0.5
# mu<- -20.25
# sig<-1
# inds<-c(0,seq_len(nx/2-1),seq(-nx/2,-1))
# xx<-inds*dx
# plot(Re(fftw::FFT(dx*dnorm(xx,mu,sig))),type='l')
# plot(Re(fftw::IFFT(exp(-2*(pi*sig*inds/(dx*nx))^2-2*pi*mu*1i*inds/(dx*nx))))~xx,type='l')
#
#
# xx<-seq(-nx*dx/2,nx*dx/2-dx,length.out=nx)
# K<-pi/(nx*dx^2)
# tmp<-exp(-2*K^2*10^2*(xx[nx/2+1]-xx[1]-abs(xx[nx/2+1]-xx))^2-2i*K*(xx-xx[1])*(0-xx[1])-log(dx))
# plot(dnorm(xx,0,10)~xx,type='l')
# plot(Re(fftw::IFFT(tmp))~xx,type='l')

#it seems the characteristic function of a BM convolution kernel doesn't quite follow the normal distribution rules since reference to
#the original scale of x is unimportant other than dx

# #more direct
# plot(Im(loc.base)*(0-xx[1]))
# #this is the expression!
# plot(Im(-1i*(xx-xx[1])*2*pi/dx^2/1024)*(0-xx[1]))
# Im(-2i*pi/1024*(xx-xx[1])*(mu-xx[1])/dx)
#
# plot(tmp)
# plot((scale.base)[1:(nx/2+1)])
# plot(-2*((xx[1:(1+nx/2)]-xx[1])*pi/(dx^2*nx))^2)
# plot((-2*(xx-xx[1])^2*pi^2/(dx^4*nx^2))[1:(nx/2+1)])
# sqrt(((xx[1:(1+nx/2)]-xx[1])*pi/dx/)^2/tmp)
#
# plot(c((-2*(xx-xx[1])^2*pi^2/(dx^4*nx^2))[1:(nx/2+1)],
#        (-2*(xx[nx]+dx-xx)^2*pi^2/(dx^4*nx^2))[(nx/2+2):nx]))
# plot(scale.base)
#
# (xx[nx]-xx[1])*(nx)/(nx-1)/2
#
#
#
# (xx[nx/2+1]-xx[1])*49/100
#
# (nx*dx^2)/(xx[nx]-xx[1])
#
# nx*dx^2
#
#
# (x*mu-x*x0-mu*x0+x0^2)/dx^2
# dx*512*dx
# diff(range(xx))/2
#
# tail(
# cbind(Im(loc.base),
#       Im(-1i*(xx-xx[1])*2*pi))
# )
# Im(loc.base)/Im(-1i*(xx-xx[1])*2*pi)
#
# tail(seq(0,2*pi/dx,length.out=nx))
# plot(Im(-1i*(xx-xx[1])*2*pi*1024/1025))
# #the above unfortunately break when x isn't sufficiently close to sampled point
# #even when sig2 is large, oddly enough
# #yeah, convolution doesn't seem to "rescue" dirac impulses that are "off-grid", so to speak
# #less beneficial than anticipated
# #interestingly, dirac deltas can be convolved with each other for very interesting results...
# #^the above seems to yield unpredictable, unreliable results, though
# #best to use these function after rounding locations appropriately (nearest xpt, etc.)
# #if you ever have to do something with shifted convolution, you may need some "ticker" mechanism which rounds but keeps track of the remainder as you propagate condtional likelihoods
# #aliasing requires zero padding, but zero padding requires non-FFT formation of the input vector...
# # test<-D(-1+0.1*dx,0.01)
# # plot(Re(IFFT(test)))
# # #zero padding in frequency domain?
# # test2<-c(test,rep(0,2048))
# # plot(Re(IFFT(test2)),type='l') #well...it got better, but not quite...
# # #wait
# # plot(Re(IFFT(test2))[seq(1,4096,2)],type='l') #whoa! but probs depends on how divided dx is
# # #yeah, would need to add A LOT of 0s depending on how much dx is being chopped up
