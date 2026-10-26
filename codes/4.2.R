setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/4.2-plot.RData")){
  p=2;n=7
  DmM=matrix(c(.5,.5,.5,0,.5,1,1/3-1/12*sqrt(7),.25,1/3-1/12*sqrt(7),.75,2/3+1/12*sqrt(7),.25,2/3+1/12*sqrt(7),.75),nrow=7,ncol=2,byrow=TRUE)
  D=DmM
  save(D,file="data/4.2-plot.RData")
}
load("data/4.2-plot.RData")
draw=function(){
  library(MASS)
  eqscplot(D[,1],D[,2],xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(0,1),ylim=c(0,1),pch=16,col="blue",main="minimax design")
  library(plotrix)
  for(i in 1:7)draw.circle(D[i,1],D[i,2],radius=1/6*(sqrt(7)-1))
}
pdf("../figures/4.2.pdf",width=4,height=4)
draw()
dev.off()
svg("../figures/4.2.svg",width=4,height=4)
draw()
dev.off()
