setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/4.22-plot.RData")){
  library(mined)
  library(zipcodeR)
  dat=na.omit(zip_code_db[,c(9,8,14)])
  ind=(dat[,2]<52 & dat[,2]>25)
  dat=dat[ind,]
  CAND=as.matrix(dat[,1:2])
  y=dat[,3]
  sizecode=.5+5*sqrt(y)/max(sqrt(y))
  D=SelectMinED(candidates=CAND,candlf=log(y),n=200,s=1)$points
  library(maximin)
  p=2
  D0=matrix(CAND[sample(1:dim(CAND)[1],1),],nrow=1,ncol=p)
  D2=CAND[maximin.cand(199,Xcand=CAND,Xorig=D0)$inds,]
  D2=rbind(D0,D2)
  sizecode=.25+log(y+1)/max(log(y+1))
  save(CAND,D2,D,sizecode,file="data/4.22-plot.RData")
}
load("data/4.22-plot.RData")
draw=function(){
  par(mfrow=c(1,2))
  plot(CAND,pch=16,col="pink",xlab="",ylab="",main="maximin design")
  points(D2,pch=16,col="blue")
  plot(CAND,pch=16,col="pink",xlab="",ylab="",main="minimum energy design")
  points(D,pch=16,col="blue")
}
pdf("../figures/4.22.pdf",width=12,height=4)
draw()
dev.off()
svg("../figures/4.22.svg",width=12,height=4)
draw()
dev.off()
