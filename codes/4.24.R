setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/4.24-plot.RData")){
  logf=function(para){
    l1=-40;u1=40;l2=-25;u2=10
    x1=l1+(u1-l1)*para[1]
    x2=l2+(u2-l2)*para[2]
    val=-.5*(x1^2/100+(x2+.03*x1^2-3)^2)
    return(val)
  }
  N.plot=300
  p1=seq(0,1,length.out=N.plot)
  p2=seq(0,1,length.out=N.plot)
  fc=matrix(0.0,N.plot,N.plot)
  for(i in 1:N.plot){
    for(j in 1:N.plot){
      fc[i,j]=exp(logf(c(p1[i],p2[j])))
    }
  }
  set.seed(8)
  p=2;n=20
  library(MaxPro)
  ini=MaxPro(MaxProLHD(n,p)$Design)$Design
  library(mined)
  res=mined(ini,logf,K_iter=5)
  D=res$points
  cand=res$cand
  ind=match(ini[,1],cand[,1])
  cand=cand[-ind,]
  save(fc,p1,p2,ini,cand,file="data/4.24-plot.RData")
}
load("data/4.24-plot.RData")
draw=function(){
  par(mfrow=c(1,2))
  library(fields)
  imagePlot(p1,p2,fc,xlab=expression(theta[1]),ylab=expression(theta[2]),col=cm.colors(5),main="initial MaxPro design")
  points(ini,pch=16,col="darkred")
  imagePlot(p1,p2,fc,xlab=expression(theta[1]),ylab=expression(theta[2]),col=cm.colors(5),main="MED points with annealing")
  points(ini,pch=16,col="darkred")
  points(cand,pch=16,col="blue")
}
pdf("../figures/4.24.pdf",width=8,height=4)
draw()
dev.off()
svg("../figures/4.24.svg",width=8,height=4)
draw()
dev.off()
