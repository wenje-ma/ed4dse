setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/8.4-plot.RData")){
  library(mixexp)
  D1=SLD(3,2)
  D2=SCD(3)
  save(D1,D2,file="data/8.4-plot.RData")
}
load("data/8.4-plot.RData")
library(Ternary)
draw=function(){
  par(mfrow=c(1,2))
  TernaryPlot(axis.labels=seq(0,1,by=0.1),atip=expression(x[1]),btip=expression(x[3]),ctip=expression(x[2]),alab="Fraction x1",blab="Fraction x3",clab="Fraction x2")
  TernaryPoints(D1,col=4,pch=16)
  TernaryPlot(axis.labels=seq(0,1,by=0.1),atip=expression(x[1]),btip=expression(x[3]),ctip=expression(x[2]),alab="Fraction x1",blab="Fraction x3",clab="Fraction x2")
  TernaryPoints(D2,col=4,pch=16)
}
pdf("../figures/8.4.pdf",width=8,height=4)
draw()
dev.off()
svg("../figures/8.4.svg",width=8,height=4)
draw()
dev.off()
