setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/6.6-plot.RData")){
  set.seed(3)
  A1=matrix(runif(4*2),nrow=4,ncol=2)
  B1=matrix(runif(4*2),nrow=4,ncol=2)
  C11=cbind(B1[,1],A1[,2])
  C12=cbind(A1[,1],B1[,2])
  D1=rbind(A1,C11,C12)
  set.seed(1)
  A2=matrix(runif(4*2),nrow=4,ncol=2)
  B2=matrix(runif(4*2),nrow=4,ncol=2)
  C21=cbind(B2[,1],A2[,2])
  C22=cbind(A2[,1],B2[,2])
  D2=rbind(A2,C21,C22)
  save(D1,D2,A1,C11,C12,A2,C21,C22,file="data/6.6-plot.RData")
}
load("data/6.6-plot.RData")
draw=function(){
  par(mfrow=c(1,2))
  plot(A1,xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(0,1),ylim=c(0,1),main="Sobol' Design (seed 3)")
  points(C11,col=2,pch=2)
  points(C12,col=3,pch=3)
  for(i in 1:4){
    arrows(D1[i,1],D1[i,2],D1[i+4,1],D1[i+4,2],col=4)
    arrows(D1[i,1],D1[i,2],D1[i+8,1],D1[i+8,2],col=4)
  }
  plot(A2,xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(0,1),ylim=c(0,1),main="Sobol' Design (seed 1)")
  points(C21,col=2,pch=2)
  points(C22,col=3,pch=3)
  for(i in 1:4){
    arrows(D2[i,1],D2[i,2],D2[i+4,1],D2[i+4,2],col=4)
    arrows(D2[i,1],D2[i,2],D2[i+8,1],D2[i+8,2],col=4)
  }
}
pdf("../figures/6.6.pdf",width=8,height=4)
draw()
dev.off()
svg("../figures/6.6.svg",width=8,height=4)
draw()
dev.off()
