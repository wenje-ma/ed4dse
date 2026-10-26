setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/9.2-plot.RData")){
  library(support)
  sp.sample=function(x,n,prob,tol=1e-6,iter.max=10){
    if(is.null(dim(x)))stop("x must be a matrix!")
    N=nrow(x)
    x.dist=as.matrix(dist(x))
    x.measure=c(x.dist%*%prob)
    idx=c(which.min(x.measure))
    if(prob[idx[1]]<1/N)idx[1]=which.max(prob)
    for(i in 2:n){
      measure=x.measure-apply(matrix(x.dist[,idx],ncol=(i-1)),1,sum)/i
      idx=c(idx,which.min(measure))
    }
    edist=2*sum(x.measure[idx])/n-sum(x.dist[idx,idx])/n^2
    iter=0
    while(TRUE){
      iter=iter+1
      for(i in 1:n){
        measure=x.measure-apply(x.dist[,idx[-i]],1,sum)/n
        idx[i]=which.min(measure)
      }
      edist.new=2*sum(x.measure[idx])/n-sum(x.dist[idx,idx])/n^2
      if((edist-edist.new)<tol||iter>iter.max){
        break
      } else {
        edist=edist.new
      }
    }
    return(idx)
  }
  L=matrix(c(0,0,1/2,0,0,1/2,1/2,1/2),ncol=2,byrow=TRUE)
  set.seed(1)
  N=20;p=2
  samp=sp(n=N,p=2,dist.str=c("uniform","normal"))$sp
  samp[,1]=1+2*samp[,1]
  samp[,2]=2+1/2*samp[,2]
  D=NULL
  for(i in 1:20){
    D=rbind(D,matrix(c(0,0,1/samp[i,1],0,0,1/samp[i,2],1/samp[i,1],1/samp[i,2]),ncol=2,byrow=TRUE))
  }
  ind=sp.sample(D,8,prob=rep(.25/20,80))
  save(L,D,ind,file="data/9.2-plot.RData")
}
load("data/9.2-plot.RData")
draw=function(){
  par(mfrow=c(1,2))
  plot(L,,pch=16,col=4,main="Locally D-optimal Design",xlim=c(0,1),ylim=c(0,1),xlab=expression(x[1]),ylab=expression(x[2]))
  text(L+.025,"2",col=2)
  plot(D,xlim=c(0,1),ylim=c(0,1),main="Robust Design",xlab=expression(x[1]),ylab=expression(x[2]))
  points(D[ind,],pch=16,col=4)
  text(.025,.025,"2",col=2)
}
pdf("../figures/9.2.pdf",width=8,height=4)
draw()
dev.off()
svg("../figures/9.2.svg",width=8,height=4)
draw()
dev.off()
