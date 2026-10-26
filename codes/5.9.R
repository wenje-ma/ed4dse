setwd("C:/Users/18904/Github/ED4DSE/codes")
if(!dir.exists("data"))dir.create("data")
if(!dir.exists("../figures"))dir.create("../figures")
if(!file.exists("data/5.9-plot.RData")){
  source("data/lib.r")
  library(mvtnorm)
  library(randtoolbox)
  logmixture=function(x){
    v=matrix(NA,nrow=5,ncol=2)
    v[1,]=c(-10,-10)
    v[2,]=c(0,16)
    v[3,]=c(13,8)
    v[4,]=c(-9,7)
    v[5,]=c(14,-14)
    v=(v+20)/40
    sigma=array(NA,dim=c(5,2,2))
    sigma[1,,]=matrix(c(2,0.6,0.6,1),nrow=2)/40^2
    sigma[2,,]=matrix(c(2,-0.4,-0.4,2),nrow=2)/40^2
    sigma[3,,]=matrix(c(2,0.8,0.8,2),nrow=2)/40^2
    sigma[4,,]=matrix(c(3,0,0,0.5),nrow=2)/40^2
    sigma[5,,]=matrix(c(2,-0.1,-0.1,2),nrow=2)/40^2
    logf=rep(0,5)
    for(i in 1:5)logf[i]=dmvnorm(x,mean=v[i,],sigma=sigma[i,,],log=TRUE)
    logaddexp(logf)-log(5)
  }
  x1=x2=seq(0,1,length.out=101)
  x.grid=expand.grid(x1,x2)
  density=matrix(exp(apply(x.grid,1,logmixture)),101,101)
  set.seed(1)
  p=2;K=20;J=10;steps=4
  ini=sobol(K,p)
  sigma=0.2
  pmc.mn=pmc(logmixture,K,J,steps,ini,sampling="random",resampling="multinomial",sigma=sigma,sigma.adapt=FALSE,visualization=F,output="weighted samples")
  pqmc.isp=pmc(logmixture,K,J,steps,ini,sampling="qmc",resampling="sp",sigma=sigma,sigma.adapt=TRUE,visualization=F,output="weighted samples")
  pmc.center=pmc.mn$center.all;pqmc.center=pqmc.isp$center.all
  pmc.samp=pmc.mn$samp.all;pqmc.samp=pqmc.isp$samp.all
  save(x1,x2,density,K,steps,pmc.center,pqmc.center,pmc.samp,pqmc.samp,file="data/5.9-plot.RData")
}
load("data/5.9-plot.RData")
draw=function(){
  par(mfrow=c(1,2))
  contour.default(x=x1,y=x2,z=density,drawlabels=FALSE,nlevels=10,main="Population Monte Carlo",xlab=expression(x[1]),ylab=expression(x[2]))
  points(pmc.center[(K*steps+1):(K*steps+K),],col="blue",pch=16)
  points(pmc.samp,col=adjustcolor(2,0.3),pch=4)
  contour.default(x=x1,y=x2,z=density,drawlabels=FALSE,nlevels=10,main="Population Quasi-Monte Carlo",xlab=expression(x[1]),ylab=expression(x[2]))
  points(pqmc.center[(K*steps+1):(K*steps+K),],col="blue",pch=16)
  points(pqmc.samp,col=adjustcolor(2,0.3),pch=4)
}
pdf("../figures/5.9.pdf",width=8,height=4)
draw()
dev.off()
svg("../figures/5.9.svg",width=8,height=4)
draw()
dev.off()
