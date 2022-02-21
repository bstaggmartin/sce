#seems scaling messed up with these guys
.get.base<-function(nx,dx){
  base<-seq(0,pi/dx,length.out=nx+1)^2
  c(base,rev(base)[-c(1,nx+1)])
}

#' @export
gauss.DFT<-function(nx,dx){
  base<- -0.5*.get.base(nx,dx)
  out<-function(t,sig2){
    exp(t*sig2*base)
  }
}

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
#' @export
VG.DFT<-function(nx,dx){
  base<-.get.base(nx,dx)
  out<-function(t,rate,jump){
    (1+base*jump/(2*rate))^(-t*rate)
  }
  out
}

# DFT<-NIG.DFT(1024,0.1)
# seed<-rep(0,2048)
# seed[512]<-1
# seed<-FFT(seed)
# plot(Re(IFFT(seed*DFT(10,1,100))),type='l')
# sum(Re(IFFT(seed*DFT(10,10,1/10))))
