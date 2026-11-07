// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <cstdlib>
#include <cmath>
using namespace Rcpp;
using namespace arma;

class ClusterDesign {
private:
  arma::mat X;
  arma::mat centers;
  double alpha;
  int Lloyd_iter_max;
  double Lloyd_tol;
  int cen_iter_max;
  double cen_tol;

public:
  arma::uvec cluster_assign;
  arma::vec distances;

  ClusterDesign(const arma::mat &data,const arma::mat &initialCenters,double alpha,
                int Lloyd_iter_max,double Lloyd_tol,
                int cen_iter_max,double cen_tol)
    : X(data),centers(initialCenters),alpha(alpha),
      Lloyd_iter_max(Lloyd_iter_max),Lloyd_tol(Lloyd_tol),
      cen_iter_max(cen_iter_max),cen_tol(cen_tol)
  {
    cluster_assign.set_size(X.n_rows);
    distances.set_size(X.n_rows);
  }

  void assignClusters(){
    int N=X.n_rows;
    for(int j=0;j < N;j++){
      rowvec xj=X.row(j);
      vec dist=sqrt(sum(square(centers.each_row()- xj),1));
      distances(j)=min(dist);
      cluster_assign(j)=dist.index_min()+ 1;
    }
  }

  static arma::rowvec Weiszfeld(const arma::mat &clusterPoints,arma::rowvec center,
                                double alpha,int iter_max=10,double tol=1e-4){
    double obj0=datum::inf;
    if(clusterPoints.n_rows == 1){
      return clusterPoints.row(0);
    }
    for(int j=0;j < iter_max;j++){
      vec dists=sqrt(sum(square(clusterPoints.each_row()- center),1));
      double obj=accu(pow(dists,alpha));
      vec w;
      if(alpha < 2){
        w=1.0 /(pow(dists,2.0 - alpha)+ 1e-16);
      } else {
        w=pow(dists,alpha - 2.0);
      }
      w=w / sum(w);
      rowvec newCenter=(w.t()* clusterPoints).eval();
      if((obj0 - obj)/ obj0 < tol) break;
      center=newCenter;
      obj0=obj;
    }
    return center;
  }

  static arma::rowvec cqAGD(const arma::mat &clusterPoints,double alpha,
                            int iter_max=1000,double tol=1e-4){
    int n=clusterPoints.n_rows;
    int p=clusterPoints.n_cols;
    int it_num=0;
    double beta=alpha - 1;
    double diff=0.0;
    rowvec curz=mean(clusterPoints,0);
    rowvec curu=curz;
    double lambdat=1.0,lambdatp1=1.0;
    rowvec prevu(p);
    double gammat=(1.0 - lambdat)/ lambdatp1;
    while(it_num < iter_max){
      lambdat=lambdatp1;
      lambdatp1=(1.0 + sqrt(1.0 + 4.0 * lambdat * lambdat))/ 2.0;
      gammat=(1.0 - lambdat)/ lambdatp1;
      prevu=curu;
      curu=curz;
      vec wts(n,fill::zeros);
      for(int i=0;i < n;i++){
        double dist=norm(curz - clusterPoints.row(i),2);
        wts(i)=(1.0 / n)*(1.0 / beta)* pow(dist,alpha - 2.0);
      }
      for(int i=0;i < n;i++){
        curu -=wts(i)*(curz - clusterPoints.row(i));
      }
      for(int j=0;j < p;j++){
        curz(j)=(1 - gammat)* curu(j)+ gammat * prevu(j);
      }
      diff=norm(curu - prevu,2);
      if(diff < tol) break;
      it_num++;
    }
    return curz;
  }

  List runLloyd(){
    int N=X.n_rows;
    assignClusters();
    double obj0=accu(pow(distances,alpha))/ N;
    double obj=obj0;
    std::vector<double> obj_hist;
    int i;
    for(i=0;i < Lloyd_iter_max;i++){
      int k=centers.n_rows;
      for(int cl=0;cl < k;cl++){
        uvec idx=find(cluster_assign ==(cl + 1));
        if(idx.n_elem > 0){
          mat clusterPoints=X.rows(idx);
          if(alpha <= 2){
            centers.row(cl)=Weiszfeld(clusterPoints,centers.row(cl),alpha,cen_iter_max,cen_tol);
          } else {
            centers.row(cl)=cqAGD(clusterPoints,alpha,cen_iter_max,cen_tol);
          }
        }
      }
      assignClusters();
      obj=accu(pow(distances,alpha))/ N;
      if((obj0 - obj)> 0 &&(obj0 - obj)/ obj0 < Lloyd_tol) break;
      obj0=obj;
      obj_hist.push_back(obj0);
    }

    Rcpp::NumericVector crit_hist_R(obj_hist.begin(),obj_hist.end());
    return List::create(
      Named("design")=centers,
      Named("cluster")=cluster_assign,
      Named("cluster_error")=obj,
      Named("total_iter")=i,
      Named("crit_hist")=crit_hist_R
);
  }
};

// [[Rcpp::export]]
double clusterError(const arma::mat& design,const arma::mat& X,double alpha=1.0){
  int N=X.n_rows;
  vec min_dist(N);

  for(int i=0;i < N;i++){
    rowvec Xi=X.row(i);
    vec distances=sqrt(sum(square(design.each_row()- Xi),1));
    min_dist(i)=min(distances);
  }

  return accu(pow(min_dist,alpha))/ N;
}

// [[Rcpp::export]]
List cluster_based_design_cpp(const arma::mat& X,arma::mat D_ini,double alpha=1.0,
                              int Lloyd_iter_max=100,double Lloyd_tol=1e-4,
                              int cen_iter_max=10,double cen_tol=1e-4){
  ClusterDesign cd(X,D_ini,alpha,Lloyd_iter_max,Lloyd_tol,cen_iter_max,cen_tol);
  return cd.runLloyd();
}

// // [[Rcpp::export]]