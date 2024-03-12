#Playing around with finite diff approx's of higher derivative using only 3 function evals
#Bottom line: doesn't seem to be possible...
xx<-seq(-5,5,length.out=100)
dx<-xx[2]-xx[1]
yy<-dnorm(xx)
inds1<-seq_along(xx)[-c(99,100)]
inds2<-seq_along(xx)[-c(1,2)]
plot(yy[-c(1,100)]~xx[-c(1,100)],type="l")
plot((yy[inds2]-yy[inds1])/(2*dx)~xx[-c(1,100)],type="l")
plot((yy[inds2]+yy[inds1]-2*yy[-c(1,100)])/dx^2~xx[-c(1,100)],type="l")
plot((yy[inds2]-yy[inds1])/(64*dx^3)~xx[-c(1,100)],type="l") #terrible approximation...
plot(diff((yy[inds2]+yy[inds1]-2*yy[-c(1,100)])/dx^2)/dx~xx[-c(1,2,100)],type="l")
plot(diff(diff((yy[inds2]+yy[inds1]-2*yy[-c(1,100)])/dx^2)/dx)/dx~xx[-c(1,2,3,100)],type="l")
plot(diff(diff(diff((yy[inds2]+yy[inds1]-2*yy[-c(1,100)])/dx^2)/dx)/dx)/dx~xx[-c(1,2,3,4,100)],type="l")

#gaussian derivatives have a pleasing form:
plot(8*(-4*xx^5+20*xx^3-15*xx)*exp(-xx^2),type="l")
plot(4*(4*xx^4-12*xx^2+3)*exp(-xx^2),type="l")
lines(4*(-2*xx^3+3*xx)*exp(-xx^2))
lines(2*(2*xx^2-1)*exp(-xx^2))
lines(2*xx*exp(-xx^2))


#Playing around with fourier transforms of step functions and the like
#(Modifying potentials to enfore reflexive rather than periodic boundary conditions)
nx<-1024
dx<-0.1
xx<-c(seq(0,dx*nx,dx),seq(-dx*(nx-1),-dx,dx))
perm.inds<-c((nx+2):(2*nx),1:(nx+1))
base<-.get.base(nx,dx)
base2<-base[[2]]
base<-base[[1]]

tmp<-xx^-3
tmp[1]<-0
plot(Im(fftw::IFFT(tmp)[perm.inds]),type="l")
tmp<-xx^-2
tmp[1]<-0
plot(Re(fftw::IFFT(tmp)[perm.inds]),type="l")
tmp<-xx^-1
tmp[1]<-0
plot(Im(fftw::IFFT(tmp)[perm.inds]),type="l")
plot(Re(fftw::IFFT(xx^0)[perm.inds]),type="l")
plot(Im(fftw::IFFT(xx)[perm.inds]),type="l")
plot(Re(fftw::FFT(xx^2)[perm.inds]),type="l")
plot(Im(fftw::FFT(xx^3)[perm.inds]),type="l")
plot(Re(fftw::FFT(xx^4)[perm.inds]),type="l")

plot(Im(fftw::FFT(-sign(xx))[perm.inds]),type="l")
plot(sign(xx))
plot(xx)

tmp<-xx^-1
tmp[1]<-1
tmp<-(sign(xx)+1)/2
plot(Im(fftw::FFT(tmp)[perm.inds]),type="l")
plot(Im(fftw::IFFT(xx)[perm.inds]),type="l")
#interesting...it's like a weirdly cut-off derivative of the dirac function

tmp<-xx^-1
tmp[seq(1,2*nx,2)]<-0 #alternating 0s for some reason...
plot(Im(fftw::FFT(tmp)[perm.inds]),type="l") #wow, so this actually creates the step function...neat!
#not sure if it practically helps with anything though--how to efficiently convolve the two functions?

tmp<-tmp/max(tmp)
tmp2<-Im(fftw::IFFT(xx))
tmp2<-tmp2/max(tmp2)
plot(tmp2,type="l")
lines(tmp,col="red")
#yeah, interesting connection, but not sure how to operationalize it...there's no intuitive way to modify differentiation in the correct way, is there?
