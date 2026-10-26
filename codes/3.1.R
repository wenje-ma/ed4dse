setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/3.1-plot.RData")){
  p=2
  n=7
  ind=c(9,45,31)
  D1=matrix(0,nrow=n,ncol=p)
  D2=matrix(0,nrow=n,ncol=p)
  D3=matrix(0,nrow=n,ncol=p)
  for(i in 1:3){
    set.seed(ind[i])
    D=matrix(runif(n*p),nrow=n,ncol=p)
    if(i==1)D1=D
    if(i==2)D2=D
    if(i==3)D3=D
  }
  save(D1,D2,D3,file="data/3.1-plot.RData")
}
load("data/3.1-plot.RData")
draw=function(){
  par(mfrow=c(1,3))
  plot(D1,xlim=c(0,1),ylim=c(0,1),xlab=expression(x[1]),ylab=expression(x[2]),pch=16,col="blue",main="Random Design 1")
  plot(D2,xlim=c(0,1),ylim=c(0,1),xlab=expression(x[1]),ylab=expression(x[2]),pch=16,col="blue",main="Random Design 2")
  plot(D3,xlim=c(0,1),ylim=c(0,1),xlab=expression(x[1]),ylab=expression(x[2]),pch=16,col="blue",main="Random Design 3")
  par(mfrow=c(1,1))
}
pdf("../figures/3.1.pdf",width=12,height=4)
draw()
dev.off()
svg("../figures/3.1.svg",width=12,height=4)
draw()
dev.off()
