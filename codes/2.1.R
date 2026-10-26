setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/2.1-plot.RData")){
  f=function(x)sin(10*pi*x)/(1+64*(x-.25)^2)+x^2
  test=seq(0,1,length=301)
  true=f(test)
  n=10
  D1=((1:n)-1)/(n-1)
  y1=f(D1)
  d=2*D1-1
  X=NULL
  for(i in 1:n)X=cbind(X,cos((i-1)*acos(d)))
  a=solve(X,y1)
  u=seq(0,1,length=301)
  yhat1=u;us=2*u-1
  for(i in 1:301){
    yhat1[i]=0
    for(j in 1:n)yhat1[i]=yhat1[i]+a[j]*cos((j-1)*acos(us[i]))
  }
  d=cos((2*(1:n)-1)/(2*n)*pi)
  D2=(d+1)/2
  y2=f(D2)
  X=NULL
  for(i in 1:n)X=cbind(X,cos((i-1)*acos(d)))
  a=solve(X,y2)
  yhat2=u;us=2*u-1
  for(i in 1:301){
    yhat2[i]=0
    for(j in 1:n)yhat2[i]=yhat2[i]+a[j]*cos((j-1)*acos(us[i]))
  }
  save(D1,y1,D2,y2,test,true,yhat1,yhat2,file="data/2.1-plot.RData")
}
load("data/2.1-plot.RData")
draw=function(){
  par(mfrow=c(1,2))
  plot(test,true,type="l",xlab="x",ylab="y",lty=2,col=1,ylim=c(min(true)-.25,max(true)+.25),main="Polynomial Interpolation (Equi-spaced)")
  points(cbind(D1,y1),pch=16,col="blue")
  lines(test,yhat1,col=3)
  legend("bottomright",legend=c("truth","prediction"),lty=c(2,1),col=c(1,3),bty="n")
  plot(test,true,type="l",xlab="x",ylab="y",lty=2,col=1,ylim=c(min(true)-.25,max(true)+.25),main="Polynomial Interpolation (Chebyshev)")
  points(cbind(D2,y2),pch=16,col="blue")
  lines(test,yhat2,col=3)
  legend("bottomright",legend=c("truth","prediction"),lty=c(2,1),col=c(1,3),bty="n")
  par(mfrow=c(1,1))
}
pdf("../figures/2.1.pdf",width=8,height=4)
draw()
dev.off()
svg("../figures/2.1.svg",width=8,height=4)
draw()
dev.off()
