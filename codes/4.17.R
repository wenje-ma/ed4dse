setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/4.17-plot.RData")){
  p=2;n=9
  set.seed(1)
  D0=cbind(rep(c(1,2,3),rep(3,3)),1:9,sample(1:9))
  D0[,2:3]=(D0[,2:3]-.5)/9
  library(SLHD)
  a=maximinSLHD(3,3,2)
  D1=a$StandDesign
  library(MaxPro)
  set.seed(3)
  ini=cbind(MaxPro(MaxProLHD(n,p)$Design)$Design,rep(1:3,c(3,3,3)))
  a=MaxProQQ(InitialDesign=ini,p_nom=1)
  a$measure
  D=a$Design
  save(D0,D1,D,p,n,file="data/4.17-plot.RData")
}
load("data/4.17-plot.RData")
draw=function(){
  par(mfrow=c(1,3))
  plot(D0[,2:3],bty="n",pch=D0[,1],col=D0[,1]+1,xlab=expression(x[1]),ylab=expression(x[2]),main="random LHD",xlim=c(0,1),ylim=c(0,1))
  for(i in 1:n){
    segments(i/n,0,i/n,1,lty=3)
    segments(0,i/n,1,i/n,lty=3)
  }
  for(i in 0:(n/3)){
    segments(3*i/n,0,3*i/n,1)
    segments(0,3*i/n,1,3*i/n)
  }
  plot(D1[,2:3],bty="n",pch=D1[,1],col=D1[,1]+1,xlab=expression(x[1]),ylab=expression(x[2]),main="SLHD",xlim=c(0,1),ylim=c(0,1))
  for(i in 1:n){
    segments(i/n,0,i/n,1,lty=3)
    segments(0,i/n,1,i/n,lty=3)
  }
  for(i in 0:(n/3)){
    segments(3*i/n,0,3*i/n,1)
    segments(0,3*i/n,1,3*i/n)
  }
  plot(D[,1:2],bty="n",pch=D[,3],col=D[,3]+1,xlab=expression(x[1]),ylab=expression(x[2]),main="MaxPro")
  for(i in 0:n){
    segments(i/n,0,i/n,1,lty=3)
    segments(0,i/n,1,i/n,lty=3)
  }
  for(i in 0:(n/3)){
    segments(3*i/n,0,3*i/n,1)
    segments(0,3*i/n,1,3*i/n)
  }
}
pdf("../figures/4.17.pdf",width=12,height=4)
draw()
dev.off()
svg("../figures/4.17.svg",width=12,height=4)
draw()
dev.off()
