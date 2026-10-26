setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/6.7-plot.RData")){
  library(MOFAT)
  p=2;l=4;s=1
  set.seed(s)
  D1=mofat(p,l,method="uniform")
  A1=D1[1:4,]
  C11=D1[(l+1):(2*l),]
  C12=D1[(2*l+1):(3*l),]
  set.seed(s)
  D2=mofat(p,l,method="projection")
  A2=D2[1:4,]
  C21=D2[(l+1):(2*l),]
  C22=D2[(2*l+1):(3*l),]
  save(D1,D2,A1,C11,C12,A2,C21,C22,file="data/6.7-plot.RData")
}
load("data/6.7-plot.RData")
draw=function(){
  par(mfrow=c(1,2))
  plot(A1,xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(0,1),ylim=c(0,1),main="MOFAT Design")
  points(C11,col=2,pch=2)
  points(C12,col=3,pch=3)
  for(i in 1:4){
    arrows(A1[i,1],A1[i,2],C11[i,1],C11[i,2],col=4,lty=1)
    arrows(A1[i,1],A1[i,2],C12[i,1],C12[i,2],col=4,lty=1)
  }
  plot(A2,xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(0,1),ylim=c(0,1),main="MOFAT Design with projections")
  points(C21,col=2,pch=2)
  points(C22,col=3,pch=3)
  for(i in 1:4){
    arrows(A2[i,1],A2[i,2],C21[i,1],C21[i,2],col=4,lty=1)
    arrows(A2[i,1],A2[i,2],C22[i,1],C22[i,2],col=4,lty=1)
  }
}
pdf("../figures/6.7.pdf",width=8,height=4)
draw()
dev.off()
svg("../figures/6.7.svg",width=8,height=4)
draw()
dev.off()
