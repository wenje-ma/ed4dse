setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/4.3-plot.RData")){
  p=2;n=7
  D1=matrix(c(0,.1,.1,0,0,1,1,1,1,.5,1,0,.5,.5),nrow=n,ncol=p,byrow=TRUE)
  D2=matrix(c(0,.1,.1,0,0,1,.9,1,1,.9,1,0,.5,.5),nrow=n,ncol=p,byrow=TRUE)
  save(D1,D2,file="data/4.3-plot.RData")
}
load("data/4.3-plot.RData")
draw=function(){
  par(mfrow=c(1,2))
  plot(D1,xlim=c(0,1),ylim=c(0,1),xlab=expression(x[1]),ylab=expression(x[2]),pch=16,col="blue",main="Index = 1")
  segments(0,.1,.1,0,col=2)
  plot(D2,xlim=c(0,1),ylim=c(0,1),xlab=expression(x[1]),ylab=expression(x[2]),pch=16,col="blue",main="Index = 2")
  segments(0,.1,.1,0,col=2)
  segments(.9,1,1,.9,col=2)
  par(mfrow=c(1,1))
}
pdf("../figures/4.3.pdf",width=8,height=4)
draw()
dev.off()
svg("../figures/4.3.svg",width=8,height=4)
draw()
dev.off()
