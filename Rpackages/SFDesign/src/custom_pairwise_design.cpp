// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <string>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

class LHDPairDesignOptimizer_custom {
protected:
  arma::mat X;
  arma::vec d;
  int n;
  int p;
  int num_passes;
  int max_iter;
  int total_iter;
  double temp;
  double decay;
  int no_update_iter_max;

  std::string optimizationMethod;

  int computePosition(int row,int h,int n){
    return(int)(row + 1 - pow((double)(h + 1),2)* 0.5 +(n - 0.5)*(h + 1)- n - 1);
  }

public:
  LHDPairDesignOptimizer_custom(const arma::mat& design,int num_passes,int max_iter,
                      double temp,double decay=0.95,int no_update_iter_max=400,
                      const std::string &method = "deterministic")
    : X(design),num_passes(num_passes),max_iter(max_iter),total_iter(0),
      temp(temp),decay(decay),
      optimizationMethod(method){
    n=design.n_rows;
    p=design.n_cols;
    this->no_update_iter_max=std::min(no_update_iter_max,5 * n *(n - 1)* p);
  }
  virtual ~LHDPairDesignOptimizer_custom(){}

  void setOptimizationMethod(const std::string &method){
    optimizationMethod=method;
  }

  virtual arma::vec computeDistanceMatrix(const arma::mat &A)=0;
  virtual double computeCriterion(const arma::vec &d)=0;
  virtual arma::vec updateDistanceMatrix(arma::mat &A,int col,int selrow1,int selrow2,arma::vec d)=0;

  List optimizeDet(){
    d=computeDistanceMatrix(X);
    double critbest=computeCriterion(d);
    std::vector<double> xcrit_hist;
    for(int pass=0;pass < num_passes;pass++){
      bool changed=false;
      for(int row1=0;row1 < n - 1;row1++){
        for(int row2=row1 + 1;row2 < n;row2++){
          for(int col=0;col < p;col++){
            total_iter++;
            if(total_iter > max_iter) goto endloop;
            arma::mat X_try=X;
            std::swap(X_try(row1,col),X_try(row2,col));
            arma::vec d_try=updateDistanceMatrix(X_try,col,row1,row2,d);
            double crit_try=computeCriterion(d_try);

            if(crit_try < critbest){
              X(row1,col)=X_try(row1,col);
              X(row2,col)=X_try(row2,col);
              critbest=crit_try;
              d=d_try;
              changed=true;
            }
            xcrit_hist.push_back(critbest);
          }
        }
      }
      if(!changed) break;
    }
    endloop:
      Rcpp::NumericVector crit_hist_R(xcrit_hist.begin(),xcrit_hist.end());
      return List::create(Named("design")=X,
                          Named("total_iter")=total_iter,
                          Named("criterion")=critbest,
                          Named("crit_hist")=crit_hist_R);
  }

  List optimizeSA(){
    d=computeDistanceMatrix(X);
    arma::mat X_best=X;
    double critbest=computeCriterion(d);
    double xcrit=critbest;
    arma::mat X_try=X;
    int ipert=1;
    bool ichange=true;
    std::vector<double> xcrit_hist;
    while(ichange){
      ichange=false;
      ipert=1;
      while(ipert < no_update_iter_max){
        if(total_iter > max_iter) break;
        total_iter++;
        int col=arma::randi(distr_param(0,p-1));
        int row1=arma::randi(distr_param(0,n-1));
        int row2=arma::randi(distr_param(0,n-2));
        if(row2 >= row1){
          row2 ++;
        }
        X_try=X;
        std::swap(X_try(row1,col),X_try(row2,col));

        arma::vec d_try=updateDistanceMatrix(X_try,col,row1,row2,d);
        double crit_try=computeCriterion(d_try);

        if(crit_try < critbest){
          ichange=true;
          X_best=X_try;
          critbest=crit_try;
          X(row1,col)=X_try(row1,col);
          X(row2,col)=X_try(row2,col);
          xcrit=crit_try;
          d=d_try;
          ipert=1;
        } else {
          ipert=ipert + 1;
          if(crit_try < xcrit){
            X(row1,col)=X_try(row1,col);
            X(row2,col)=X_try(row2,col);
            d=d_try;
            xcrit=crit_try;
            ichange=true;
          } else if(arma::randu()< exp(-(crit_try - xcrit)/ temp)){
            X(row1,col)=X_try(row1,col);
            X(row2,col)=X_try(row2,col);
            d=d_try;
            xcrit=crit_try;
          }
        }
        xcrit_hist.push_back(xcrit);
      }
      temp=temp * decay;
    }

    Rcpp::NumericVector xcrit_hist_R(xcrit_hist.begin(),xcrit_hist.end());
    return List::create(Named("design")=X_best,
                        Named("total_iter")=total_iter,
                        Named("crit_hist")=xcrit_hist_R,
                        Named("criterion")=critbest);
  }

  List optimize(){
    if(optimizationMethod == "deterministic"){
      return optimizeDet();
    } else if(optimizationMethod == "sa"){
      return optimizeSA();
    } else {
      Rcpp::Rcout << "Unknown optimization method: " << optimizationMethod
                  << ". Using deterministic." << std::endl;
      return optimizeDet();
    }
  }
};

class CustomLHDOptimizer : public LHDPairDesignOptimizer_custom {
private:
  std::function<arma::vec(const arma::mat&)> user_computeDistanceMatrix;
  std::function<double(const arma::vec&)> user_computeCriterion;
  std::function<arma::vec(arma::mat&,int,int,int,arma::vec)> user_updateDistanceMatrix;

public:
  CustomLHDOptimizer(Rcpp::Function r_computeDistanceMatrix,
                    Rcpp::Function r_computeCriterion,
                    Rcpp::Function r_updateDistanceMatrix,
                    const arma::mat &design,int num_passes,int max_iter,
                    double temp=0,double decay=0.95,int no_update_iter_max=400,
                    const std::string &method = "deterministic")
    : LHDPairDesignOptimizer_custom(design,num_passes,max_iter,
      temp,decay,no_update_iter_max,method)
  {
    user_computeDistanceMatrix=[r_computeDistanceMatrix](const arma::mat &A)-> arma::vec {
      Rcpp::NumericMatrix A_rcpp=Rcpp::wrap(A);
      Rcpp::NumericVector result=r_computeDistanceMatrix(A_rcpp);
      return Rcpp::as<arma::vec>(result);
    };
    user_computeCriterion=[r_computeCriterion](const arma::vec &d)-> double {
      Rcpp::NumericVector d_rcpp=Rcpp::wrap(d);
      Rcpp::NumericVector result=r_computeCriterion(d_rcpp);
      return result[0];
    };
    user_updateDistanceMatrix=[r_updateDistanceMatrix](arma::mat &A,int col,int selrow1,int selrow2,arma::vec d)-> arma::vec {
      Rcpp::NumericMatrix A_rcpp=Rcpp::wrap(A);
      Rcpp::NumericVector d_rcpp=Rcpp::wrap(d);
      Rcpp::NumericVector result=r_updateDistanceMatrix(A_rcpp,col,selrow1,selrow2,d_rcpp);
      return Rcpp::as<arma::vec>(result);
    };
  }
  arma::vec computeDistanceMatrix(const arma::mat &A) override {
    return user_computeDistanceMatrix(A);
  }
  double computeCriterion(const arma::vec &d) override {
    return user_computeCriterion(d);
  }
  arma::vec updateDistanceMatrix(arma::mat &A,int col,int selrow1,int selrow2,arma::vec d) override {
    return user_updateDistanceMatrix(A,col,selrow1,selrow2,d);
  }
};

// [[Rcpp::export]]
List customLHDOptimizer_cpp(Rcpp::Function r_computeDistanceMatrix,
                            Rcpp::Function r_computeCriterion,
                            Rcpp::Function r_updateDistanceMatrix,
                            arma::mat design,int num_passes=10,int max_iter=1e6,
                            double temp=0,double decay=0.95,int no_update_iter_max=400,
                            std::string method = "deterministic"){
  CustomLHDOptimizer optimizer(r_computeDistanceMatrix,
                                r_computeCriterion,
                                r_updateDistanceMatrix,
                                design,num_passes,max_iter,
                                temp,decay,no_update_iter_max,
                                method);
  return optimizer.optimize();
}