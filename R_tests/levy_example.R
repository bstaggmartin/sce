tree<-pbtree(n=50,scale=1)
dat<-sim.levy(tree,model=1,rate=0.5,jump=2^2,sig2=0.1^2)
phenogram(tree,dat[,1],ftype='off')
levy<-make.levy(tree,dat[,1],model=1,nx=4096)
levy[['lik.fun']](log(c(0.5,2^2,0.1^2))) #oof, infinite likelihood under true params...had to bump up to nx=2048
tmp<-optim(rnorm(3),levy[['lik.fun']],control=list(trace=TRUE),method='BFGS')
#definitely need to check DFTs against simulations to ensure proper scaling, but it actually looks pretty good
#the likelihoods seem to be on very different scales across levy models, unsure what that's about
test<-levy[['recon.fun']](tmp$par)
xpts<-attr(test,'xpts')
dx<-diff(xpts[1:2])
anc.means<-dx*.rowSums(sweep(test,2,xpts,'*'),dim(test)[1],dim(test)[2])
nodes<-get('nodes',envir=environment(levy[['recon.fun']]),inherits=FALSE)[,2]
ord<-order(nodes)
anc.means<-anc.means[ord]
names(anc.means)<-c(tree[['tip.label']],seq_len(tree[['Nnode']])+Ntip(tree))
rbind(anc.means[tree$tip.label],dat[,1]) #looks right
phenogram(tree,anc.means,ftype='off')
tips<-get('tips',envir=environment(levy[['recon.fun']]),inherits=FALSE)
matplot(t(test[-tips,]),type='l',x=xpts) #some bimodality for compound poisson! cool!
