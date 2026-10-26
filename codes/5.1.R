setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/5.1-plot.RData")){
  set.seed(1)
  n=7
  x1=sort(runif(n))
  d1=cbind(x1,(1:n)/n)
  x2=((1:n)-.5)/n
  d2=cbind(x2,(1:n)/n)
  save(d1,d2,x1,x2,n,file="data/5.1-plot.RData")
}
load("data/5.1-plot.RData")
draw=function(){
  par(mfrow=c(1,2))
  plot(d1,xlim=c(0,1),ylim=c(0,1),pch=16,col="blue",xlab="x",ylab="Distribution Function",main="Random")
  segments(0,0,x1[1],0,col=2)
  for(i in 1:(n-1))segments(x1[i],i/n,x1[i+1],i/n,col=2)
  segments(x1[n],1,1,1,col=2)
  abline(0,1,col=3,lty=2)
  plot(d2,xlim=c(0,1),ylim=c(0,1),pch=16,col="blue",xlab="x",ylab="Distribution Function",main="Optimal")
  segments(0,0,x2[1],0,col=2)
  for(i in 1:(n-1))segments(x2[i],i/n,x2[i+1],i/n,col=2)
  segments(x2[n],1,1,1,col=2)
  abline(0,1,col=3,lty=2)
}
pdf("../figures/5.1.pdf",width=8,height=4)
draw()
dev.off()
svg("../figures/5.1.svg",width=8,height=4)
draw()
dev.off()
