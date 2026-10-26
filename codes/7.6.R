setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/7.6-plot.RData")){
  library(OSFD)
  library(MaxPro)
  f=function(x){
    y1=1/(x[1]^2+x[2]^2+.01)^(1/2)
    if(x[2]==0)y2=0 else if(x[1]==0)y2=pi/2 else y2=atan(x[2]/x[1])
    return(c(y1=y1,y2=y2))
  }
  set.seed(1)
  p=2
  q=2
  n0=10
  n=50
  D0=MaxProLHD(n0,p)$Design
  osfd=OSFD(D=D0,f=f,p=p,q=q,n=n,method="Greedy")
  D=osfd$D
  Y=osfd$Y
  save(D,Y,n0,n,file="data/7.6-plot.RData")
}
load("data/7.6-plot.RData")
draw=function(){
  par(mfrow=c(1,2))
  plot(D,pch=NA,main="Input Space",xlab=expression(x[1]),ylab=expression(x[2]))
  points(D[1:n0,],pch=16,col=4)
  points(D[(n0+1):n,],pch=4,col=2)
  plot(Y,pch=NA,main="Output Space",xlab=expression(y[1]),ylab=expression(y[2]))
  points(Y[1:n0,],pch=16,col=4)
  points(Y[(n0+1):n,],pch=4,col=2)
}
pdf("../figures/7.6.pdf",width=8,height=4)
draw()
dev.off()
svg("../figures/7.6.svg",width=8,height=4)
draw()
dev.off()
