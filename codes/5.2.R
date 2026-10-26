setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/5.2-plot.RData")){
  p=2;n=20
  set.seed(6)
  library(spacefillr)
  library(SFDesign)
  R=matrix(runif(n*2),ncol=2)
  S=generate_sobol_set(n,p,seed=sample(1:10000,1))
  D=uniform.optim(uniformLHD(n,p)$design)$design
  save(R,S,D,file="data/5.2-plot.RData")
}
load("data/5.2-plot.RData")
draw=function(){
  par(mfrow=c(1,3))
  plot(R,pch=16,col=4,main="Monte Carlo",xlab=expression(x[1]),ylab=expression(x[2]))
  plot(S,pch=16,col=4,main="Sobol Sequence",xlab=expression(x[1]),ylab=expression(x[2]))
  plot(D,pch=16,col=4,main="Uniform Design",xlab=expression(x[1]),ylab=expression(x[2]))
}
pdf("../figures/5.2.pdf",width=12,height=4)
draw()
dev.off()
svg("../figures/5.2.svg",width=12,height=4)
draw()
dev.off()
