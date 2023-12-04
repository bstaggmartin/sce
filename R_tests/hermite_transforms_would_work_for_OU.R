library(EQL)
xx<-seq(-10,10,length.out=100)
plot(hermite(xx,1)~xx,type="l")

x0<- -10
dx<-0.05
xpts<-seq(x0,-x0,dx)
nx<-length(xpts)
alpha<-0.2
sig2<-1
theta<-0
tmp<-seq(2*dx*alpha/sig2*(x0+dx/2-theta),by=2*dx^2*alpha/sig2,length.out=nx)
tmp<-list(exp(tmp[-nx]),exp(-tmp[-1]))
Q<-matrix(0,nx,nx)
diag(Q[-1,])<-tmp[[1]]*sig2/(2*dx^2)
diag(Q[,-1])<-tmp[[2]]*sig2/(2*dx^2)
diag(Q)<- -rowSums(Q)

#eigenvectors do seem to consist of hermite polynomials...
tmp<-eigen(Q)
invtmp<-solve(tmp$vectors)
image(tmp$vectors%*%diag(exp(1*tmp$values))%*%invtmp)
image(tmp$vectors%*%diag(exp(1*tmp$values)))



get.P<-function(xpts,alpha,theta,sig2,t){
  nx<-length(xpts)
  out<-matrix(dnorm(rep(xpts,each=nx),
                    rep(xpts,nx)*exp(-alpha*t)+theta*(1-exp(-alpha*t)),
                    sig2/(2*alpha)*(1-exp(-2*alpha*t))),
              nx,nx)
  out<-out/rowSums(out)
  out
}
n<-3
alpha<-0.5
theta<-5
sig2<-1
t<-4
P<-get.P(xpts,alpha,theta,sig2,t)
hh<-hermite((xpts-theta)*sqrt(2*alpha/sig2),n)
plot(P%*%hh,type="l",col="red")
lines(hh)
lines(hh*exp(-n*alpha*t),col="blue")
#wow...it does seem to work...
#so you can do state-dependent OU if there was a fast a way to hermite and inverse hermite transform probability density functions...
#but I'm not sure that works...
#hmmm...and it turns out Hermite transforms have an extra e^(-x^2) component that causes individual components to NOT solve OU process
#but that component seems to ensure the polynomials are orthogonal such that they can approximate arbitrary functions...
#so this is probably a dead-end, but nonetheless interesting...
#there also seems to be some necessary phase-shifting when theta is non-0, which could be problematic...
#of course, xpts could always be transformed to ensure theta is at 0 for both discrete states?

#yeah, found a transform that seems to do the trick with varying parameters
#(still get some differences between numerical solution and theoretical solution, but it could just be weird boundary conditions...)
#In any case, it still doesn't help all that much because you need to break down pdfs in terms of hermite polynomials, which doesn't seem to be straight-forward at all

hh<-Reduce("+",lapply(0:6,function(ii) (-1)^(ii+1)*exp(-ii)*hermite((xpts-theta)*sqrt(2*alpha/sig2),ii)))
plot(hh~xpts,type="l",ylim=c(-100,100))
#yeah, seems really hard to make sensible distributions out of these functions...

#looking into it, it looks like this all still might work with the hermite transform because you can simply multiply the original function by x^2 as well...
#but then it wouldn't be following a strict OU process, right?
#but some papers seem to indicate that you can use this trick intermediately to still break down a function into a sum of hermite polynomials

#would seem to require working with OU in a modified space or something...
