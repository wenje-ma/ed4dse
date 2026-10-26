setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/4.8-plot.RData")){
  p=2;n=7
  library(SFDesign)
  D=maximinLHD(n,p)$design
  save(D,n,file="data/4.8-plot.RData")
}
load("data/4.8-plot.RData")
draw=function(){
  plot(D,pch=16,axes=FALSE,bty="n",col="blue",xlim=c(0,1),ylim=c(0,1),xlab=expression(x[1]),ylab=expression(x[2]),main="Maximin LHD",asp=1)
  axis(1,at=seq(0,1,by=1/7),labels=round(seq(0,1,by=1/7),1))
  axis(2,at=seq(0,1,by=1/7),labels=round(seq(0,1,by=1/7),1))
  for(i in 1:(n+1)){
    segments(0,(i-1)/n,1,(i-1)/n)
    segments((i-1)/n,0,(i-1)/n,1)
  }
}
pdf("../figures/4.8.pdf",width=4,height=4)
draw()
dev.off()
svg("../figures/4.8.svg",width=4,height=4)
draw()
dev.off()
