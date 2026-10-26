setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/4.14-plot.RData")){
  branin=function(x){
    x1=x[1]*15-5
    x2=x[2]*15
    val1=(x2-5.1/(4*pi^2)*(x1^2)+5/pi*x1-6)^2
    val2=10*(1-1/(8*pi))*cos(x1)+10
    return(val1+val2)
  }
  p=2;n=10
  set.seed(123)
  library(MaxPro)
  D=MaxPro(MaxProLHD(n,p)$Design)$Design
  CAND=CandPoints(N=10000,2)
  A=MaxProAugment(ExistDesign=D,CandDesign=CAND,nNew=10)$Design[(n+1):20,]
  y=apply(D,1,branin)
  library(rkriging)
  a=Fit.Kriging(D,y,kernel.parameters=list(type="Gaussian"))
  true=apply(A,1,branin)
  pred=Predict.Kriging(a,A)$mean
  rmse0=sqrt(mean((pred-true)^2))
  rmse=numeric(100)
  for(i in 1:100){
    u=matrix(runif(p*10),nrow=10,ncol=p)
    pred=Predict.Kriging(a,u)$mean
    true=apply(u,1,branin)
    rmse[i]=sqrt(mean((pred-true)^2))
  }
  N.plot=250
  p1=seq(0,1,length=N.plot)
  p2=seq(0,1,length=N.plot)
  fc=matrix(0,N.plot,N.plot)
  for(i in 1:N.plot){
    for(j in 1:N.plot)fc[i,j]=branin(c(p1[i],p2[j]))
  }
  save(fc,p1,p2,rmse,D,A,rmse0,file="data/4.14-plot.RData")
}
load("data/4.14-plot.RData")
draw=function(){
  par(mar=c(5,5,4,2))
  par(mfrow=c(1,2))
  library(fields)
  imagePlot(p1,p2,fc,xlab=expression(x[1]),ylab=expression(x[2]),col=cm.colors(12,rev=TRUE),main="Validation points")
  points(D,pch=16,col="blue")
  points(A,pch=2,col="darkred")
  boxplot(rmse,main="RMSE")
  abline(h=rmse0,col="darkred")
  text(.8,rmse0+2,"sequential MaxPro")
}
pdf("../figures/4.14.pdf",width=8,height=4)
draw()
dev.off()
svg("../figures/4.14.svg",width=8,height=4)
draw()
dev.off()
