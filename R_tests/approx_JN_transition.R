#the idea here is to avoid having to integrate over all possible numbers of jumps
#instead, just divide into no jump vs yes jump categories
#calculate avg number of jumps in yes jump category to approximate yes jump dist as a single normal rather than mixture of normals
#seems to work particularly well for small timesteps, unsurprisingly
#so, you could solve JN on sphere by simply solving the heat equation twice!
#once for no jump dist, then yes jump dist, then combine with appropriate probability
approx.JN<-function(sigma2,mu,lambda,tau2,delta){
  function(x,x0,t){
    noj<-exp(-lambda*t)
    cnoj<-1-noj
    avgj<-lambda*t/cnoj
    noj*dnorm(x,x0+t*mu,sqrt(t*sigma2))+cnoj*dnorm(x,x0+avgj*delta+t*mu,sqrt(avgj*tau2+t*sigma2))
  }
}
#really, that's not too bad! Particularly over short time intervals...
true.JN<-function(sigma2,mu,lambda,tau2,delta){
  function(x,x0,t){
    max.j<-qpois(1-1e-16,lambda*t)
    njs<-0:max.j
    pjs<-dpois(njs,lambda*t)
    colSums(
      matrix(
        pjs*
          dnorm(rep(x,each=max.j+1),
                x0+t*mu+njs*delta,
                sqrt(t*sigma2+njs*tau2)),
        max.j+1,
        length(x)
        )
      )
  }
}
sigma2<-1
mu<- -1
JN.fun1<-approx.JN(sigma2,mu,0.1,100,2)
JN.fun2<-true.JN(sigma2,mu,0.1,100,2)
xx<-seq(-20,20,length.out=1000)
t<-10
tmp1<-dnorm(xx,0+t*mu,sqrt(t*sigma2),log=TRUE)
tmp2<-log(JN.fun2(xx,0,t))
tmp3<-log(JN.fun1(xx,0,t))
plot(tmp1~xx,type="l",ylim=c(min(tmp2,tmp3),max(tmp1)))
lines(tmp2~xx,col="blue")
lines(tmp3~xx,col="red")
