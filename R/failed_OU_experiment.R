x0<- -10
dx<-0.1
xpts<-seq(x0,-x0,dx)
nx<-length(xpts)
alpha<-0.01
sig2<-0.1
theta<- -5
# tmp<-alpha/sig2*xpts^2-2*alpha*theta/sig2*xpts
#alpha/sig2*((xpts+dx)^2-2*theta*(xpts+dx)-xpts^2+2*theta*xpts)
#xpts^2+2*xpts*dx+dx^2-2*theta*xpts-2*theta*dx-xpts^2+2*theta*xpts
#2*xpts*dx+dx^2-2*theta*dx
#2*dx*alpha/sig2*(xpts+dx/2-theta)
tmp<-seq(2*dx*alpha/sig2*(x0+dx/2-theta),by=2*dx^2*alpha/sig2,length.out=nx)
tmp<-list(exp(tmp[-1]),exp(-tmp[-nx]))
Q<-matrix(0,nx,nx)
diag(Q[-1,])<-tmp[[1]]*sig2/(2*dx^2)
diag(Q[,-1])<-tmp[[2]]*sig2/(2*dx^2)
diag(Q)<- -rowSums(Q)

eig<-eigen(Q)
matplot(Re(eig$vectors),type="l")
image(Re(eig$vectors))
image(Im(eig$vectors))
plot(Re(eig$values))
plot(Im(eig$values),type="l")
plot(Re(eig$vectors[,201]),type="l")
#hmmm...do seem to be polynomials of increasing degree as you go to smaller and smaller eigenvalues
#not the same shape as the standard Hermite polynomials on wikipedia, though...
matplot(eig$vectors[,201:191],type="l",lty=1,col=rainbow(10))
#but this seems promising...

theta<-5
tmp<-seq(2*dx*alpha/sig2*(x0+dx/2-theta),by=2*dx^2*alpha/sig2,length.out=nx)
tmp<-list(exp(tmp[-nx]),exp(-tmp[-1]))
Q<-matrix(0,nx,nx)
diag(Q[-1,])<-tmp[[1]]
diag(Q[,-1])<-tmp[[2]]
diag(Q)<- -rowSums(Q)
eig2<-eigen(Q)

plot(Re(eig$values))
plot(Re(eig2$values)) #way more stable with theta=5 over theta=-5 for some reason...
#anyways...
Q<-matrix(c(-1,1,1,-1),2,2)
Q<-lapply(seq_len(nx),function(ii) {diag(Q)<-diag(Q)+c(eig$values[ii],eig2$values[ii]);Q})
tmp<-do.call(rbind,lapply(Q,function(ii) {tmp<-eigen(ii);0.5*colSums(Re(tmp$vectors%*%diag(exp(100*tmp$values))%*%solve(tmp$vectors)))}))
test1<-eig$vectors%*%diag(tmp[,1])%*%solve(eig$vectors)
plot(test1[201,]~xpts)
test2<-eig2$vectors%*%diag(tmp[,2])%*%solve(eig2$vectors)
plot(test2[1,]~xpts) #doesn't seem to work unfortunately
#worth a shot
