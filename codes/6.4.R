setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/6.4-plot.RData")){
  set.seed(10)
  par(mfrow=c(1,2))
  x=((1:4)-1)/3
  A=cbind(x,sample(x))
  delta=1/3
  C11=C12=A
  for(i in 1:4){
    C11[i,1]=A[i,1]+delta
    if(C11[i,1]<0||C11[i,1]>1)C11[i,1]=A[i,1]-delta
    C12[i,1]=C11[i,1]
    C12[i,2]=C11[i,2]+delta
    if(C12[i,2]<0||C12[i,2]>1)C12[i,2]=C11[i,2]-delta
  }
  D=rbind(A,C11,C12)
  dim(unique(D))[1]
  C21=C22=A
  for(i in 1:4){
    C21[i,1]=A[i,1]+delta
    if(C21[i,1]<0||C21[i,1]>1)C21[i,1]=A[i,1]-delta
    C22[i,2]=C21[i,2]+delta
    if(C22[i,2]<0||C22[i,2]>1)C22[i,2]=C21[i,2]-delta
  }
  D=rbind(A,C21,C22)
  dim(unique(D))[1]
  out=cbind(A,C11,C12,C21,C22)
  colnames(out)=c("A1","A2","C11.1","C11.2","C12.1","C12.2","C21.1","C21.2","C22.1","C22.2")
  save(A,C11,C12,C21,C22,file="data/6.4-plot.RData")
}
load("data/6.4-plot.RData")
draw=function(){
  par(mfrow=c(1,2))
  plot(A,xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(0,1),ylim=c(0,1),main="Morris Design (Strict OFAT)")
  points(C11,col=2,pch=2)
  points(C12,col=3,pch=3)
  for(i in 1:4){
    arrows(A[i,1],A[i,2],C11[i,1],C11[i,2],col=4,lty=1)
    arrows(C11[i,1],C11[i,2],C12[i,1],C12[i,2],col=4,lty=1)
  }
  plot(A,xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(0,1),ylim=c(0,1),main="Morris Design (Standard OFAT)")
  points(C21,col=2,pch=2)
  points(C22,col=3,pch=3)
  for(i in 1:4){
    arrows(A[i,1],A[i,2],C21[i,1],C21[i,2],col=4,lty=1)
    arrows(A[i,1],A[i,2],C22[i,1],C22[i,2],col=4,lty=1)
  }
}
pdf("../figures/6.4.pdf",width=8,height=4)
draw()
dev.off()
svg("../figures/6.4.svg",width=8,height=4)
draw()
dev.off()
