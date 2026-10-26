setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
library(MASS)
library(plotrix)
if(!file.exists("data/4.4-plot.RData")){
  p=2;n=7
  set.seed(1)
  library(SFDesign)
  D=maximinLHD(n,p)$design
  D=maximin.optim(D,sa=TRUE,find.best.ini=TRUE)$design
  r=min(dist(D))/2
  save(D,r,n,file="data/4.4-plot.RData")
}
load("data/4.4-plot.RData")
draw=function(){
  eqscplot(D[,1],D[,2],xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(0,1),ylim=c(0,1),pch=16,col="blue",main="Maximin Design")
  for(i in 1:n)draw.circle(D[i,1],D[i,2],radius=r)
}
pdf("../figures/4.4.pdf",width=4,height=4)
draw()
dev.off()
svg("../figures/4.4.svg",width=4,height=4)
draw()
dev.off()
