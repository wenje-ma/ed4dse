setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/5.6-plot.RData")){
  p=2;n=50;N.plot=300
  set.seed(1)
  library(support)
  D1=sp(n,p,dist.str=rep("normal",p))$sp
  f1=function(x)exp(-.5*sum(x^2))/(2*pi)
  p11=seq(-3,3,length.out=N.plot)
  p12=seq(-3,3,length.out=N.plot)
  fc1=matrix(apply(expand.grid(p11,p12),1,f1),N.plot,N.plot)
  dist.param=vector("list",p)
  for(l in 1:p)dist.param[[l]]=c(2,1)
  D2=sp(n,p,dist.str=rep("gamma",p),dist.param=dist.param)$sp
  f2=function(x)prod(dgamma(x,2,1))
  p21=seq(0,6,length.out=N.plot)
  p22=seq(0,6,length.out=N.plot)
  fc2=matrix(apply(expand.grid(p21,p22),1,f2),N.plot,N.plot)
  save(D1,D2,p11,p12,p21,p22,fc1,fc2,file="data/5.6-plot.RData")
}
load("data/5.6-plot.RData")
draw=function(){
  par(mfrow=c(1,2))
  image(p11,p12,fc1,xlab=expression(x[1]),ylab=expression(x[2]),col=cm.colors(5),main="Normal Distribution (Support Points)",asp=1)
  contour(p11,p12,fc1,add=TRUE,nlevels=10)
  points(D1,pch=16,col="blue")
  image(p21,p22,fc2,xlab=expression(x[1]),ylab=expression(x[2]),col=cm.colors(5),main="Gamma Distribution (Support Points)",asp=1)
  contour(p21,p22,fc2,add=TRUE,nlevels=10)
  points(D2,pch=16,col="blue")
}
pdf("../figures/5.6.pdf",width=8,height=4)
draw()
dev.off()
svg("../figures/5.6.svg",width=8,height=4)
draw()
dev.off()
