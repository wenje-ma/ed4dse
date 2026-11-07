// [[Rcpp::depends(BH)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <float.h>
#include <boost/any.hpp>
#include <time.h>
#include <sstream>
#include <string>
#include <random>

using namespace Rcpp;
using namespace std;
using namespace arma;

// // [[Rcpp::export]]
// // [[Rcpp::export]]

// [[Rcpp::export]]
double herdobj_seq(NumericVector& xx,NumericMatrix& des,NumericMatrix& distsamp,double sigma){

  double sigma2=pow(sigma,2.0);
  int pp=des.ncol();
  int nn=des.nrow();
  int NN=distsamp.nrow();
  double runsum=0.0;
  double tmp=0.0;

  for(int m=0;m<NN;m++){
    tmp=0.0;
    for(int l=0;l<pp;l++){
      tmp +=pow(xx(l)- distsamp(m,l),2.0);
    }
    runsum +=(1.0/(double) NN)* exp(-1.0/(2.0*sigma2)* tmp);
  }

  for(int i=0;i<nn;i++){
    tmp=0.0;
    for(int l=0;l<pp;l++){
      tmp +=pow(xx(l)- des(i,l),2.0);
    }
    runsum +=-(1.0/((double) nn+1.0))* exp(-1.0/(2.0*sigma2)* tmp);
  }

  return(runsum);
}

// [[Rcpp::export]]
NumericVector herdgrad_seq(NumericVector& xx,NumericMatrix& des,NumericMatrix& distsamp,double sigma){

  double sigma2=pow(sigma,2.0);
  int pp=des.ncol();
  int nn=des.nrow();
  int NN=distsamp.nrow();
  NumericVector runvec(pp,0.0);

  double tmp=0.0;

  for(int m=0;m<NN;m++){
    tmp=0.0;
    for(int l=0;l<pp;l++){
      tmp +=pow(xx(l)- distsamp(m,l),2.0);
    }
    for(int l=0;l<pp;l++){
      runvec(l)+=-(1.0/(double) NN)*(1.0/sigma2)* exp(-1.0/(2.0*sigma2)* tmp)*(xx(l)-distsamp(m,l));
    }
  }

  for(int i=0;i<nn;i++){
    tmp=0.0;
    for(int l=0;l<pp;l++){
      tmp +=pow(xx(l)- des(i,l),2.0);
    }
    for(int l=0;l<pp;l++){
      runvec(l)+=(1.0/((double) nn+1.0))*(1.0/sigma2)* exp(-1.0/(2.0*sigma2)* tmp)*(xx(l)-des(i,l));
    }
  }

  return(runvec);
}

// [[Rcpp::export]]
double herdobj_full(NumericVector& xx,NumericMatrix& des,int idx,NumericMatrix& distsamp,double sigma){

  double sigma2=pow(sigma,2.0);
  int pp=des.ncol();
  int nn=des.nrow();
  int NN=distsamp.nrow();
  double runsum=0.0;
  double tmp=0.0;

  for(int m=0;m<NN;m++){
    tmp=0.0;
    for(int l=0;l<pp;l++){
      tmp +=pow(xx(l)- distsamp(m,l),2.0);
    }
    runsum +=(1.0/(double) NN)* exp(-1.0/(2.0*sigma2)* tmp);
  }

  for(int i=0;i<nn;i++){
    if(i != idx){
      tmp=0.0;
      for(int l=0;l<pp;l++){
        tmp +=pow(xx(l)- des(i,l),2.0);
      }
      runsum +=-(1.0/((double) nn))* exp(-1.0/(2.0*sigma2)* tmp);
    }
  }

  return(runsum);
}

// [[Rcpp::export]]
NumericVector herdgrad_full(NumericVector& xx,NumericMatrix& des,int idx,NumericMatrix& distsamp,double sigma){

  double sigma2=pow(sigma,2.0);
  int pp=des.ncol();
  int nn=des.nrow();
  int NN=distsamp.nrow();
  NumericVector runvec(pp,0.0);

  double tmp=0.0;

  for(int m=0;m<NN;m++){
    tmp=0.0;
    for(int l=0;l<pp;l++){
      tmp +=pow(xx(l)- distsamp(m,l),2.0);
    }
    for(int l=0;l<pp;l++){
      runvec(l)+=-(1.0/(double) NN)*(1.0/sigma2)* exp(-1.0/(2.0*sigma2)* tmp)*(xx(l)-distsamp(m,l));
    }
  }

  for(int i=0;i<nn;i++){
    if(i != idx){
      tmp=0.0;
      for(int l=0;l<pp;l++){
        tmp +=pow(xx(l)- des(i,l),2.0);
      }
      for(int l=0;l<pp;l++){
        runvec(l)+=1.0/((double) nn)*(1.0/sigma2)* exp(-1.0/(2.0*sigma2)* tmp)*(xx(l)-des(i,l));
      }
    }
  }

  return(runvec);
}

// [[Rcpp::export]]
double pspobj_seq(NumericVector& xx,NumericMatrix& des,NumericMatrix& distsamp,
                    double lambda,double nu){

  int pp=des.ncol();
  int nn=des.nrow();
  int NN=distsamp.nrow();
  double runsum=0.0;
  double tmp=0.0;

  for(int m=0;m<NN;m++){
    tmp=1.0;
    for(int l=0;l<pp;l++){
      tmp=tmp *(1.0 /(pow(xx(l)- distsamp(m,l),2.0)+ lambda));
    }

    tmp=pow(tmp,nu)*pow(lambda,nu*pp);
    runsum +=(1.0/(double) NN)* tmp;
  }

  for(int i=0;i<nn;i++){
    tmp=1.0;
    for(int l=0;l<pp;l++){
      tmp=tmp *(1.0 /(pow(xx(l)- des(i,l),2.0)+ lambda));
    }
    tmp=pow(tmp,nu)*pow(lambda,nu*pp);
    runsum +=-(1.0/((double) nn+1.0))* tmp;
  }

  return(runsum);
}

// [[Rcpp::export]]
NumericVector pspgrad_seq(NumericVector& xx,NumericMatrix& des,NumericMatrix& distsamp,
                            double lambda,double nu){

  int pp=des.ncol();
  int nn=des.nrow();
  int NN=distsamp.nrow();
  NumericVector runvec(pp,0.0);

  double tmp=1.0;
  for(int m=0;m<NN;m++){
    tmp=1.0;
    for(int l=0;l<pp;l++){
      tmp=tmp *(1.0 /(pow(xx(l)- distsamp(m,l),2.0)+ lambda));
    }
    tmp=(1.0/(double) NN)* pow(tmp,nu)*pow(lambda,nu*pp);
    for(int l=0;l<pp;l++){
      runvec(l)+=tmp /(pow(xx(l)- distsamp(m,l),2.0)+ lambda)*(xx(l)- distsamp(m,l))*(- 2.0* nu);
    }
  }

  for(int i=0;i<nn;i++){
    tmp=1.0;
    for(int l=0;l<pp;l++){
      tmp=tmp *(1.0 /(pow(xx(l)- des(i,l),2.0)+ lambda));
    }
    tmp=-(1.0/((double) nn+1.0))* pow(tmp,nu)*pow(lambda,nu*pp);
    for(int l=0;l<pp;l++){
      runvec(l)+=tmp /(pow(xx(l)- des(i,l),2.0)+ lambda)*(xx(l)- des(i,l))*(- 2.0*nu);
    }
  }

  return(runvec);
}

// [[Rcpp::export]]
double pspobj_full(NumericVector& xx,NumericMatrix& des,int idx,NumericMatrix& distsamp,
                  double lambda,double nu){

  int pp=des.ncol();
  int nn=des.nrow();
  int NN=distsamp.nrow();
  double runsum=0.0;
  double tmp=0.0;

  for(int m=0;m<NN;m++){
    tmp=1.0;
    for(int l=0;l<pp;l++){
      tmp=tmp *(1.0 /(pow(xx(l)- distsamp(m,l),2.0)+ lambda));
    }

    tmp=pow(tmp,nu)*pow(lambda,nu*pp);
    runsum +=(1.0/(double) NN)* tmp;
  }

  for(int i=0;i<nn;i++){
    tmp=1.0;
    if(i != idx){
      for(int l=0;l<pp;l++){
        tmp=tmp *(1.0 /(pow(xx(l)- des(i,l),2.0)+ lambda));
      }
      tmp=pow(tmp,nu)*pow(lambda,nu*pp);
      runsum +=-(1.0/(double) nn)* tmp;
    }
  }

  return(runsum);
}

// [[Rcpp::export]]
NumericVector pspgrad_full(NumericVector& xx,NumericMatrix& des,int idx,NumericMatrix& distsamp,
                          double lambda,double nu){

  int pp=des.ncol();
  int nn=des.nrow();
  int NN=distsamp.nrow();
  NumericVector runvec(pp,0.0);

  double tmp=1.0;
  for(int m=0;m<NN;m++){
    tmp=1.0;
    for(int l=0;l<pp;l++){
      tmp=tmp *(1.0 /(pow(xx(l)- distsamp(m,l),2.0)+ lambda));
    }
    tmp=(1.0/(double) NN)* pow(tmp,nu)*pow(lambda,nu*pp);
    for(int l=0;l<pp;l++){
      runvec(l)+=tmp /(pow(xx(l)- distsamp(m,l),2.0)+ lambda)*(xx(l)- distsamp(m,l))*(- 2.0* nu);
    }
  }

  for(int i=0;i<nn;i++){
    tmp=1.0;
    if(i != idx){
      for(int l=0;l<pp;l++){
        tmp=tmp *(1.0 /(pow(xx(l)- des(i,l),2.0)+ lambda));
      }
      tmp=-(1.0/(double) nn)* pow(tmp,nu)*pow(lambda,nu*pp);
      for(int l=0;l<pp;l++){
        runvec(l)+=tmp /(pow(xx(l)- des(i,l),2.0)+ lambda)*(xx(l)- des(i,l))*(- 2.0*nu);
      }
    }
  }

  return(runvec);
}

// [[Rcpp::export]]
double pspobj_seq2(NumericVector& xx,NumericMatrix& des,NumericMatrix& distsamp,
                  double lambda,double nu){

  int pp=des.ncol();
  int nn=des.nrow();
  int NN=distsamp.nrow();
  double runsum=0.0;
  double tmp=0.0;

  for(int m=0;m<NN;m++){
    tmp=1.0;
    for(int l=0;l<pp;l++){
      tmp=tmp *(1.0 /(abs(xx(l)- distsamp(m,l))+ lambda));
    }

    tmp=pow(tmp,nu)*pow(lambda,nu*pp);
    runsum +=(1.0/(double) NN)* tmp;
  }

  for(int i=0;i<nn;i++){
    tmp=1.0;
    for(int l=0;l<pp;l++){
      tmp=tmp *(1.0 /(abs(xx(l)- des(i,l))+ lambda));
    }
    tmp=pow(tmp,nu)*pow(lambda,nu*pp);
    runsum +=-(1.0/((double) nn+1.0))* tmp;
  }

  return(runsum);
}

// [[Rcpp::export]]
NumericVector pspgrad_seq2(NumericVector& xx,NumericMatrix& des,NumericMatrix& distsamp,
                          double lambda,double nu){

  int pp=des.ncol();
  int nn=des.nrow();
  int NN=distsamp.nrow();
  NumericVector runvec(pp,0.0);

  double tmp=1.0;
  for(int m=0;m<NN;m++){
    tmp=1.0;
    for(int l=0;l<pp;l++){
      tmp=tmp *(1.0 /(abs(xx(l)- distsamp(m,l))+ lambda));
    }
    tmp=(1.0/(double) NN)* pow(tmp,nu)*pow(lambda,nu*pp);
    for(int l=0;l<pp;l++){
      runvec(l)+=tmp /(abs(xx(l)- distsamp(m,l))+ lambda)*(xx(l)- distsamp(m,l))/ abs(xx(l)- distsamp(m,l))*(- nu);
    }
  }

  for(int i=0;i<nn;i++){
    tmp=1.0;
    for(int l=0;l<pp;l++){
      tmp=tmp *(1.0 /(abs(xx(l)- des(i,l))+ lambda));
    }
    tmp=-(1.0/((double) nn+1.0))* pow(tmp,nu)*pow(lambda,nu*pp);
    for(int l=0;l<pp;l++){
      runvec(l)+=tmp /(abs(xx(l)- des(i,l))+ lambda)*(xx(l)- des(i,l))/ abs(xx(l)- des(i,l))*(- nu);
    }
  }

  return(runvec);
}

// [[Rcpp::export]]
double pspobj_full2(NumericVector& xx,NumericMatrix& des,int idx,NumericMatrix& distsamp,
                    double lambda,double nu){

  int pp=des.ncol();
  int nn=des.nrow();
  int NN=distsamp.nrow();
  double runsum=0.0;
  double tmp=0.0;

  for(int m=0;m<NN;m++){
    tmp=1.0;
    for(int l=0;l<pp;l++){
      tmp=tmp *(1.0 /(abs(xx(l)- distsamp(m,l))+ lambda));
    }

    tmp=pow(tmp,nu)*pow(lambda,nu*pp);
    runsum +=(1.0/(double) NN)* tmp;
  }

  for(int i=0;i<nn;i++){
    tmp=1.0;
    if(i != idx){
      for(int l=0;l<pp;l++){
        tmp=tmp *(1.0 /(abs(xx(l)- des(i,l))+ lambda));
      }
      tmp=pow(tmp,nu)*pow(lambda,nu*pp);
      runsum +=-(1.0/(double) nn)* tmp;
    }
  }

  return(runsum);
}

// [[Rcpp::export]]
NumericVector pspgrad_full2(NumericVector& xx,NumericMatrix& des,int idx,NumericMatrix& distsamp,
                            double lambda,double nu){

  int pp=des.ncol();
  int nn=des.nrow();
  int NN=distsamp.nrow();
  NumericVector runvec(pp,0.0);

  double tmp=1.0;
  for(int m=0;m<NN;m++){
    tmp=1.0;
    for(int l=0;l<pp;l++){
      tmp=tmp *(1.0 /(abs(xx(l)- distsamp(m,l))+ lambda));
    }
    tmp=(1.0/(double) NN)* pow(tmp,nu)*pow(lambda,nu*pp);
    for(int l=0;l<pp;l++){
      runvec(l)+=tmp /(abs(xx(l)- distsamp(m,l))+ lambda)*(xx(l)- distsamp(m,l))/ abs(xx(l)- distsamp(m,l))*(- nu);
    }
  }

  for(int i=0;i<nn;i++){
    tmp=1.0;
    if(i != idx){
      for(int l=0;l<pp;l++){
        tmp=tmp *(1.0 /(abs(xx(l)- des(i,l))+ lambda));
      }
      tmp=-(1.0/(double) nn)* pow(tmp,nu)*pow(lambda,nu*pp);
      for(int l=0;l<pp;l++){
        runvec(l)+=tmp /(abs(xx(l)- des(i,l))+ lambda)*(xx(l)- des(i,l))/ abs(xx(l)- des(i,l))*(- nu);
      }
    }
  }

  return(runvec);
}