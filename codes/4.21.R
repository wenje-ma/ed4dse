setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/4.21-plot.RData")){
  n1=100;n2=10;p=2
  set.seed(1)
  library(SFDesign)
  D1=maxpro.optim(maxproLHD(n1,p)$design)$design
  maxpro.crit(D1)
  D2=maxpro.remove(D1,n.remove=n1-n2,delta=1/n1^2)
  maxpro.crit(D2)
  save(D1,D2,file="data/4.21-plot.RData")
}
load("data/4.21-plot.RData")
draw=function(){
  plot(D1,xlim=c(0,1),ylim=c(0,1),pch=16,col="blue",main="MaxPro-based Nested Design",xlab=expression(x[1]),ylab=expression(x[2]))
  points(D2,pch=2,col=2)
  legend("topright",legend=c("Full (n=100)","Subset (n=10)"),pch=c(16,2),col=c("blue",2),bty="n")
}
pdf("../figures/4.21.pdf",width=4,height=4)
draw()
dev.off()
svg("../figures/4.21.svg",width=4,height=4)
draw()
dev.off()
