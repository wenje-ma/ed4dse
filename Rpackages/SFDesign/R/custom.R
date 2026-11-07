customLHD=function(compute.distance.matrix,compute.criterion,update.distance.matrix,n,p,design=NULL,max.sa.iter=1e6,temp=0,decay=0.95,no.update.iter.max=400,num.passes=10,max.det.iter=1e6,method="full",scaled=TRUE){
  if(is.null(design)){
    design=randomLHD(n,p)
  }
  if(!scaled){
    design=(apply(design,2,rank)-0.5)/n
  }
  if(method=="deterministic"){
    result=customLHDOptimizer_cpp(compute.distance.matrix,compute.criterion,update.distance.matrix,design,num.passes,max.det.iter,temp,decay,no.update.iter.max,method)
  }else if(method=="sa"){
    result=customLHDOptimizer_cpp(compute.distance.matrix,compute.criterion,update.distance.matrix,design,num.passes,max.sa.iter,temp,decay,no.update.iter.max,method)
  }else if(method=="full"){
    result=customLHDOptimizer_cpp(compute.distance.matrix,compute.criterion,update.distance.matrix,design,num.passes,max.sa.iter,temp,decay,no.update.iter.max,"sa")
    crit_hist=result$crit_hist
    total_iter=result$total_iter
    result=customLHDOptimizer_cpp(compute.distance.matrix,compute.criterion,update.distance.matrix,result$design,num.passes,max.det.iter,temp,decay,no.update.iter.max,"deterministic")
    result$crit_hist=c(crit_hist,result$crit_hist)
    result$total_iter=c(total_iter,result$total_iter)
  }
  return(result)
}