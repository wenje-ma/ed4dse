setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/9.3-plot.RData")){
  best=read.table(file="data/txts/run51.txt")$V1
  pdet=c(seq(0,5000,length=2001),seq(5000,120000,length.out=4601)[-1])
  sorp=read.csv(file="data/exp_sorption_130c_8.7torr_483nm.csv")
  desorp=read.csv(file="data/exp_desorption_130c_8.7torr_483nm.csv")
  correct=35358.89275-34644.0896-1
  yexp=c(sorp[,2],desorp[,2]+correct)
  texp=c(sorp[,1]^2,desorp[,1]^2+62973)
  save(best,pdet,yexp,texp,file="data/9.3-plot.RData")
} else {
  load("data/9.3-plot.RData")
}
draw=function(){
  plot(sqrt(texp),yexp,main="Model Prediction vs Real Data",xlab=expression(sqrt(t)),ylab="mass uptake",type="l",lty=2)
  lines(sqrt(pdet),best,col=2)
  legend("bottomright",legend=c("measured","predicted"),col=c(1,2),lty=c(2,1),bty="n")
}
pdf("../figures/9.3.pdf",width=4,height=4)
draw()
dev.off()
svg("../figures/9.3.svg",width=4,height=4)
draw()
dev.off()
