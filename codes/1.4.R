setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/1.4-plot.RData")){
  x=seq(0,1,length=301)
  n=15
  D1=((1:n)-.5)/n
  y1=sin(10*pi*x)/(1+64*(x-.5)^2)
  D2=qbeta(((1:n)-.5)/n,.5,.5)
  y2=1/(1+64*(x-.5)^2)
  n1=n*2/3; n2=n/3
  D0=((1:n2)-.5)/n2
  D3=c(D0,qbeta(((1:n1)-.5)/n1,2,10))
  y3=2/3*dbeta(x,2,10)+1/3*dbeta(x,10,2)
  set.seed(1)
  x0=rbeta(1000*n1,2,10)
  x0=c(x0,rbeta(1000*n2,10,2))
  library(support)
  D4=sp(n,1,dist.samp=cbind(x0))$sp
  y4=2/3*dbeta(x,2,10)+1/3*dbeta(x,10,2)
  save(D1,D2,D3,D4,x,y1,y2,y3,y4,file="data/1.4-plot.RData")
}
load("data/1.4-plot.RData")
draw=function(){
  par(mar=c(1,1,1,1))
  par(mfrow=c(2,2))
  plot(x,y1,type="l",axes=FALSE,xlab="",ylab="",ylim=c(min(y1)-.05,max(y1)),main="Approximate a Rough Function")
  points(cbind(D1,min(y1)-.05),pch=16,col="blue")
  plot(x,y2,type="l",axes=FALSE,xlab="",ylab="",ylim=c(min(y2)-.05,max(y2)),main="Approximate a Smooth Function")
  points(cbind(D2,min(y2)-.05),pch=16,col="blue")
  plot(x,y3,type="l",axes=FALSE,xlab="",ylab="",ylim=c(min(y3)-.05,max(y3)),main="Optimization")
  points(cbind(D3,min(y3)-.05),pch=16,col="blue")
  plot(x,y4,type="l",axes=FALSE,xlab="",ylab="",ylim=c(min(y4)-.05,max(y4)),main="Uncertainty Propagation")
  points(cbind(D4,min(y4)-.05),pch=16,col="blue")
}
pdf("../figures/1.4.pdf",width=8,height=8)
draw()
dev.off()
svg("../figures/1.4.svg",width=8,height=8)
draw()
dev.off()
