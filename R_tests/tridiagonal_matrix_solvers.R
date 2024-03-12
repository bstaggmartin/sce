solve.tridiag.periodic<-function(a,b,c,d,n){
  gamma<-5 #random choice
  b[1]<-b[1]-gamma
  b[n]<-b[n]-c[n]*a[1]/gamma
  u<-c(gamma,rep(0,n-2),c[n])
  v<-c(1,rep(0,n-2),a[1]/gamma)
  y<-solve.tridiag(a[-1],b,c[-n],d,n)
  q<-solve.tridiag(a[-1],b,c[-n],u,n)
  y-q*crossprod(v,y)[1,1]/(1+crossprod(v,q)[1,1])
}

solve.tridiag<-function(a,b,c,d,n){
  out<-numeric(n)
  for(i in seq_len(n-1)){
    ii<-i+1
    w<-a[i]/b[i]
    b[ii]<-b[ii]-w*c[i]
    d[ii]<-d[ii]-w*d[i]
  }
  out[n]<-d[n]/b[n]
  for(i in (n-1):1){
    out[i]<-(d[i]-c[i]*out[i+1])/b[i]
  }
  out
}
test<-diag(rnorm(100,0,0.1))
diag(test[,-1])<-rnorm(99,0,1)
diag(test[-1,])<-rnorm(99,0,1)
soln<-rnorm(100,0,1)
solve(test,soln)
solve.tridiag(diag(test[-1,]),diag(test),diag(test[,-1]),soln,100) #seems pretty stable...maybe as n gets large?
#well, might be "stable enough" for your purpose in any case!

a<-rnorm(100)
b<-rnorm(100,0,0.1)
c<-rnorm(100)
d<-rnorm(100)
test<-solve.tridiag.periodic(a,b,c,d,100)
mat<-diag(b)
diag(mat[-1,])<-a[-1]
mat[1,n]<-a[1]
diag(mat[,-1])<-c[-n]
mat[n,1]<-c[n]
plot(test)
plot(solve(mat,d))
mat%*%solve(mat,d)
mat%*%test
