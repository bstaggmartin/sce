
#3/8: variance gamma process simulation and fitting now concordant!
t<-1
freq<-0
jump<-1
skew<-0
yy<-rgamma(1e6,t*freq,freq)
test<-rnorm(1e6,skew*yy,sqrt(jump*yy))
tmp.hist<-hist(test,freq=FALSE,breaks=100)

xx<-seq(-6,6,length.out=1024)
dx<-xx[2]-xx[1]
tmp.res<-Re(fftw::IFFT(get.DFTs(1024,dx,"VG")[[1]](t,freq,jump,skew)))[c(1538:2048,1:513)]
tmp.res<-tmp.res/max(tmp.res)*max(tmp.hist$density)
plot(tmp.res~xx) #hmmm...something is going wrong here


.rinvgauss<-function(n,mean,shape){
  nu<-rnorm(n)
  nu.sq<-nu^2
  mean.sq<-mean^2
  meannu.sq<-mean.sq*nu.sq
  shape2<-shape*2
  x<-mean+meannu.sq/shape2-mean/shape2*sqrt(4*mean*shape*nu.sq+meannu.sq*nu.sq)
  tmp.inds<-!is.nan(x)
  nn<-sum(tmp.inds)
  inds<-vector(length=n)
  tmp<-rep(mean,length.out=n)[tmp.inds]
  inds[tmp.inds]<-runif(nn)>tmp/(tmp+x[tmp.inds])
  x[inds]<-rep(mean.sq,length.out=n)[inds]/x[inds]
  x[!tmp.inds]<-rep(shape,length.out=n)[!tmp.inds]/rchisq(n-nn,1)
  x
}



#there's some problem here...but I'm having trouble figuring it out
#not perfect but very close...
#there seems to be something wrong with the characteristic function, but I'm not sure what...
#something to do with the translation with the bases...multiplying by sqrt(freq) seems to kinda rescue things...
#something wrong with skew still maybe?...but we're getting closer
#FINALLY figured it out based on Schraiber's pulsr code...
#gah, what a waste of half a week
#I wonder what the hell is up with all the different wrong-looking parameterizations online...

#3/11: below is DONE!
#to do: better parameterization with less interparam dependence
#i.e., skew from -1 to 1, make jump a function of freq
#also make .rinvgauss reduce to inverse chisq when mean is Inf!
#but then NIG should be good to go
#and you get the cauchy process to boot!

t<-1
freq<-0
jump<-0.2
skew<- 0.95

# nn<-1e6
# hist(rnorm(nn,0,sqrt(.rinvgauss(nn,1,sqrt(freq^2-skew^2))))+
#        rnorm(nn,0,sqrt(.rinvgauss(nn,1,sqrt(freq^2-skew^2)))))
# hist(rnorm(nn,0,sqrt(.rinvgauss(nn,2,2*sqrt(freq^2-skew^2)))))

# hist(sqrt(jump*.rinvgauss(1e6,1/freq,jump*sqrt(1-skew^2))))
# hist(sqrt(jump/freq*.rinvgauss(1e6,1,jump*freq*sqrt(1-skew^2))))

# yy<-.rinvgauss(1e6,1,jump*freq*sqrt(1-skew^2))
# test<-rnorm(1e6,
#             yy*sqrt(freq)*skew*jump^2,
#             sqrt(jump/freq*yy))
# yy<-.rinvgauss(1e6,1/freq,jump*sqrt(1-skew^2))
# test<-rnorm(1e6,
#             yy*sign(skew)*freq*sqrt(abs(skew))*jump,
#             sqrt(jump*yy))

# yy<-if(freq==0){
#   #reduces to cauchy random variable!
#   1/rgamma(10,0.5,jump^2*t^2/2)
#   #OR
#   # jump^2*t^2/rgamma(10,0.5,0.5)
#   #OR
#   # jump^2*t^2/rchisq(10,1)
# }else{
#   .rinvgauss(1e6,jump*t/sqrt(freq^2-skew^2),jump^2*t^2)
# }

yy<-.rinvgauss(1e6,jump*t/sqrt(1-skew^2),freq^2*jump^2*t^2)
test<-rnorm(1e6,
            yy*freq*skew,
            sqrt(yy))
plot(density(test,from=-10,to=10))

xx<-seq(-10,10,length.out=1024)
dx<-xx[2]-xx[1]
lines(Re(fftw::IFFT(get.DFTs(1024,dx,"NIG")[[1]](t,freq,jump,skew)))[c(1538:2048,1:513)]/dx~xx,
     type="l")

plot(Re(fftw::IFFT(get.DFTs(1024,dx,"NIG")[[1]](t,1,0.5,-0.9)))[c(1538:2048,1:513)]/dx~xx,
      type="l")
