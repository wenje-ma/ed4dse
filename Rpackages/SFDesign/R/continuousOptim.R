continuous.optim=function(D.ini,objective,gradient=NULL,iteration=10,sa=FALSE,sa.objective=NULL){
  D0=as.matrix(D.ini)
  n=nrow(D0)
  p=ncol(D0)
  objective.list.flag=is.list(objective(D.ini))
  if(is.null(gradient)&(!objective.list.flag)){
    optim.obj=objective
    D1=D0
    for(i in 1:iteration){
      a=stats::optim(D1,optim.obj,method='L-BFGS-B',lower=rep(0,n*p),upper=rep(1,n*p))
      D1=matrix(a$par,nrow=n,ncol=p)
    }
  }else{
    if(objective.list.flag){
      optim.obj=objective
    }else{
      optim.obj=function(x){
        return(list("objective"=objective(x),"gradient"=gradient(x)))
      }
    }
    D1=D0
    for(i in 1:iteration){
      a=nloptr::nloptr(c(D1),optim.obj,opts=list("algorithm"="NLOPT_LD_LBFGS","maxeval"=100),lb=rep(0,n*p),ub=rep(1,n*p))
      D1=matrix(a$sol,nrow=n,ncol=p)
    }
  }
  if(sa){
    if(is.null(sa.objective)){
      warning("No sa objective is provided!")
    }
    result=GenSA::GenSA(D1,function(x) sa.objective(x),lower=rep(0,n*p),upper=rep(1,n*p))
    D1=matrix(result$par,ncol=p)
  }
  return(D1)
}