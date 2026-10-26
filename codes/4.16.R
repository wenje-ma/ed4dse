setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/4.16-plot.RData")){
  constraint=function(x){
    c1=(x[1]-sqrt(50*(x[2]-0.52)^2+2)+1)
    c2=(sqrt(120*(x[2]-0.48)^2+1)-0.75-x[1])
    c3=(0.65^2-x[1]^2-x[2]^2)
    return(c(c1,c2,c3))
  }
  x1=x2=matrix(NA,nrow=3,ncol=1001)
  x2.seq=seq(0,1,length.out=1001)
  x2[1,]=x2.seq
  x1[1,]=sqrt(50*(x2.seq-0.52)^2+2)-1
  x2[2,]=x2.seq
  x1[2,]=sqrt(120*(x2.seq-0.48)^2+1)-0.75
  x2[3,]=x2.seq
  x1[3,]=sqrt(0.65^2-x2.seq^2)
  N=1e5
  set.seed(1)
  library(lhs)
  lhs.all=randomLHS(N,2)
  lhs.gval=t(apply(lhs.all,1,constraint))
  lhs.out.idx=apply(lhs.gval,1,function(x)return(any(x>0)))
  lhs.feasible=lhs.all[!lhs.out.idx,]
  n=20
  library(MaxPro)
  a=MaxProAugment(ExistDesign=lhs.feasible[1,],CandDesign=lhs.feasible[-1,],nNew=n-1)
  D=a$Design
  all=lhs.all;feasible=lhs.feasible;maxpro=D
  save(x1,x2,all,feasible,lhs.feasible,maxpro,file="data/4.16-plot.RData")
}
load("data/4.16-plot.RData")
draw=function(){
  par(mfrow=c(1,3))
  plot(all,col="green",pch=18,xlim=c(0,1),ylim=c(0,1),xlab=expression(x[1]),ylab=expression(x[2]),main="initial points")
  for(i in 1:3)lines(x1[i,],x2[i,])
  plot(feasible,col="red",pch=16,xlim=c(.3,.8),ylim=c(.35,.55),xlab=expression(x[1]),ylab=expression(x[2]),main="candidate points")
  for(i in 1:3)lines(x1[i,],x2[i,])
  plot(maxpro,col="blue",pch=16,xlim=c(.3,.8),ylim=c(.35,.55),xlab=expression(x[1]),ylab=expression(x[2]),main="MaxPro design")
  for(i in 1:3)lines(x1[i,],x2[i,])
}
pdf("../figures/4.16.pdf",width=12,height=4)
draw()
dev.off()
svg("../figures/4.16.svg",width=12,height=4)
draw()
dev.off()
