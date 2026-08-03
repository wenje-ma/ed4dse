// [[Rcpp::depends(BH)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <float.h>
#include <boost/any.hpp>
#include <time.h>
#include <limits.h>
#include <sstream>
#include <string>
#include <random>

using namespace Rcpp;
using namespace std;
using namespace arma;

double m_eps=pow(10.0,-10.0);

void flush_console(){
#if !defined(WIN32)&& !defined(__WIN32)&& !defined(__WIN32__)
  R_FlushConsole();
#endif
}

// [[Rcpp::export]]
void printBar(double prop){

  int barWidth=70;

  flush_console();
  Rcout << "[";
  int pos=barWidth * prop;
  for(int i=0;i < barWidth;++i){
    if(i < pos) Rcout << "=";
    else if(i == pos) Rcout << ">";
    else Rcout << " ";
  }
  Rcout << "] " << int(prop * 100.0)<< " %\r";

}

namespace sftrabbit {

template <typename RealType=double>
class beta_distribution
{
public:
  typedef RealType result_type;

  class param_type
  {
  public:
    typedef beta_distribution distribution_type;

    explicit param_type(RealType a=2.0,RealType b=2.0)
      : a_param(a),b_param(b){ }

    RealType a() const { return a_param;}
    RealType b() const { return b_param;}

    bool operator==(const param_type& other) const
    {
      return(a_param == other.a_param &&
              b_param == other.b_param);
    }

    bool operator!=(const param_type& other) const
    {
      return !(*this == other);
    }

  private:
    RealType a_param,b_param;
  };

  explicit beta_distribution(RealType a=2.0,RealType b=2.0)
    : a_gamma(a),b_gamma(b){ }
  explicit beta_distribution(const param_type& param)
    : a_gamma(param.a()),b_gamma(param.b()){ }

  void reset(){ }

  param_type param() const
  {
    return param_type(a(),b());
  }

  void param(const param_type& param)
  {
    a_gamma=gamma_dist_type(param.a());
    b_gamma=gamma_dist_type(param.b());
  }

  template <typename URNG>
  result_type operator()(URNG& engine)
  {
    return generate(engine,a_gamma,b_gamma);
  }

  template <typename URNG>
  result_type operator()(URNG& engine,const param_type& param)
  {
    gamma_dist_type a_param_gamma(param.a()),
    b_param_gamma(param.b());
    return generate(engine,a_param_gamma,b_param_gamma);
  }

  result_type min() const { return 0.0;}
  result_type max() const { return 1.0;}

  result_type a() const { return a_gamma.alpha();}
  result_type b() const { return b_gamma.alpha();}

  bool operator==(const beta_distribution<result_type>& other) const
  {
    return(param()== other.param()&&
            a_gamma == other.a_gamma &&
            b_gamma == other.b_gamma);
  }

  bool operator!=(const beta_distribution<result_type>& other) const
  {
    return !(*this == other);
  }

private:
  typedef std::gamma_distribution<result_type> gamma_dist_type;

  gamma_dist_type a_gamma,b_gamma;

  template <typename URNG>
  result_type generate(URNG& engine,
                        gamma_dist_type& x_gamma,
                        gamma_dist_type& y_gamma)
  {
    result_type x=x_gamma(engine);
    return x /(x + y_gamma(engine));
  }
};

template <typename CharT,typename RealType>
std::basic_ostream<CharT>& operator<<(std::basic_ostream<CharT>& os,
                                      const beta_distribution<RealType>& beta)
{
  os << "~Beta(" << beta.a()<< "," << beta.b()<< ")";
  return os;
}

template <typename CharT,typename RealType>
std::basic_istream<CharT>& operator>>(std::basic_istream<CharT>& is,
                                      beta_distribution<RealType>& beta)
{
  std::string str;
  RealType a,b;
  if(std::getline(is,str,'(')&& str == "~Beta" &&
      is >> a && is.get()== ',' && is >> b && is.get()== ')'){
    beta=beta_distribution<RealType>(a,b);
  } else {
    is.setstate(std::ios::failbit);
  }
  return is;
}

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector csample_num(NumericVector x,
                            int size,
                            bool replace){
  NumericVector prob(x.size());
  for(int i=0;i<x.size();i++){
    prob[i]=1.0/x.size();
  }
  NumericVector ret=Rcpp::RcppArmadillo::sample(x,size,replace,prob);
  return(ret);
}

double sgn(double val){
  int ret=(val > 0.0)-(val < 0.0);
  return((double) ret);
}

// [[Rcpp::export]]
double energycrit(NumericMatrix& Rcpp_point,NumericMatrix& Rcpp_des){

  double ret=0.0;

  int dim_num=Rcpp_point.ncol();
  int point_num=Rcpp_point.nrow();
  int des_num=Rcpp_des.nrow();

    double runcrit=0.0;
    for(int i=0;i<des_num;i++){
      for(int j=0;j<point_num;j++){
        double runcrittmp=0.0;
        for(int k=0;k<dim_num;k++){
          runcrittmp +=pow(Rcpp_point(j,k)- Rcpp_des(i,k),2.0);
        }
        runcrit +=sqrt(runcrittmp);
      }
    }
    runcrit=runcrit *(2.0/(des_num*point_num));

    double runcrit2=0.0;
    for(int i=0;i<des_num;i++){
      for(int j=0;j<des_num;j++){
        double runcrittmp2=0.0;
        for(int k=0;k<dim_num;k++){
          runcrittmp2 +=pow(Rcpp_des(i,k)- Rcpp_des(j,k),2.0);
        }
        runcrit2 +=sqrt(runcrittmp2);
      }
    }
    runcrit2=runcrit2 /((double)(des_num*des_num));

  ret=runcrit - runcrit2;
  return(ret);

}

// [[Rcpp::export]]
NumericVector starL2cpp(NumericMatrix& D,NumericMatrix& cn,int num_proc){

  int n=D.nrow();
  int dimen=cn.ncol();
  NumericVector ret(cn.nrow());

  for(int l=0;l<cn.nrow();l++){
    double term1=pow(3.0,-dimen);
    double runsum=0.0;
    for(int i=0;i<n;i++){
      double runprod=1.0;
      for(int j=0;j<dimen;j++){
        runprod=runprod *(1 - pow(D(i,cn(l,j)-1),2));
      }
      runsum +=runprod;
    }
    double term2=(pow(2.0,1-dimen)/n)* runsum;

    runsum=0.0;
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        double runprod=1.0;
        for(int k=0;k<dimen;k++){
          runprod=runprod *(1-max(D(i,cn(l,k)-1),D(j,cn(l,k)-1)));
        }
        runsum +=runprod;
      }
    }
    double term3=(1/pow(n,2.0))*runsum;

    ret(l)=sqrt(term1-term2+term3);
  }

  return(ret);

}

// [[Rcpp::export]]
NumericVector energy_norm_cpp(NumericMatrix& yMat,NumericMatrix& cn,int num_proc)
{

  int n=yMat.nrow();
  int d=cn.ncol();
  NumericVector retvec(cn.nrow());

  for(int l=0;l<cn.nrow();l++){
    std::vector<double> y(d*n);
    for(int i=0;i<n;i++){
      for(int j=0;j<d;j++){
        y[i*d+j]=yMat(i,cn(l,j)-1);
      }
    }
    int    i,j,k,p,maxterms=2000;
    double D=(double) d;
    double meanyy,meanzz;
    double delta,eps=1.0e-7;
    double normy,yy,dif,sum,sum0,term;
    double lg0,lg1,logak,loggk;

    lg0=lgamma(D/2.0);
    lg1=lgamma((D+1.0)/2.0);
    meanzz=2.0 * exp(lg1 - lg0);

    arma::vec meanyz(n);
    for(i=0;i<n;i++){
      yy=0.0;
      for(j=0;j<d;j++){
        dif=y[i*d+j]* y[i*d+j];
        yy +=dif;
      }
      normy=sqrt(yy);
      delta=1.0;
      sum=0.0;
      k=0;
      while(delta > eps && k < maxterms){
        sum0=sum;
        logak=(k+1)*log(yy)- lgamma(k+1)- k*M_LN2 -
          log(2*k+1)- log(2*k+2);
        loggk=lg1 + lgamma(k+1.5)- lgamma(k+D/2+1);
        term=exp(logak + loggk);
        if(k % 2 == 0)
          sum +=term;
        else
          sum -=term;
        delta=fabs(sum - sum0);
        k++;
      }
      if(delta < eps)
        meanyz(i)=meanzz/M_SQRT2 + M_SQRT_2dPI * sum;
      else {
        meanyz(i)=normy;
      }
    }

    double ret=0.0;
    arma::vec normyvec(n);
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        yy=0.0;
        for(k=0;k<d;k++){
          dif=pow(y[i*d+k]- y[j*d+k],2.0);
          yy +=dif;
        }
        normy=sqrt(yy);
        ret +=meanyz(i)+ meanyz(j)- meanzz - normy;
      }
    }
    retvec(l)=ret / pow((double) n,2.0);
  }

  return(retvec);
}

// [[Rcpp::export]]
double obj_qsp(arma::vec& des,arma::mat& distsamp,double q){
  int NN=distsamp.n_rows;
  int pp=distsamp.n_cols;
  int nn=des.n_elem/pp;
  double runsum=0.0;
  double runsum2=0.0;
  double tmp=0.0;
  double ret=0.0;

  for(int i=0;i<nn;i++){
    for(int m=0;m<NN;m++){
      tmp=0.0;
      for(int j=0;j<pp;j++){
        tmp +=pow(des(i*pp+j)-distsamp(m,j),2);
      }
      runsum +=pow(tmp,q/2.0);
    }
  }
  runsum=2.0/((double) nn * NN)* runsum;

  for(int i=0;i<nn;i++){
    for(int k=0;k<nn;k++){
      tmp=0.0;
      for(int j=0;j<pp;j++){
        tmp +=pow(des(i*pp+j)-des(k*pp+j),2);
      }
      runsum2 +=pow(tmp,q/2.0);
    }
  }
  runsum2=runsum2/((double) pow(nn,2.0));

  ret=runsum - runsum2;
  return(ret);
}

// [[Rcpp::export]]
arma::vec grad_qsp(arma::vec& des,arma::mat& distsamp,double q){
  int NN=distsamp.n_rows;
  int pp=distsamp.n_cols;
  int nn=des.n_elem/pp;
  double tmp=0.0;
  arma::vec ret(nn*pp);
  arma::vec ret2(nn*pp);
  for(int i=0;i<(nn*pp);i++){
    ret(i)=0.0;
    ret2(i)=0.0;
  }
  arma::vec tmpvec(pp);

  for(int i=0;i<nn;i++){
    for(int m=0;m<NN;m++){

      for(int j=0;j<pp;j++){
        tmpvec(j)=0.0;
      }
      tmp=0.0;

      for(int j=0;j<pp;j++){
        tmpvec(j)+=des(i*pp+j)-distsamp(m,j);
        tmp +=pow(tmpvec(j),2);
      }
      tmpvec=tmpvec / pow(tmp,(2.0-q)/2.0);

      for(int j=0;j<pp;j++){
        ret(i*pp+j)+=tmpvec(j);
      }

    }
  }
  ret=(2.0 * q)/((double) nn * NN)* ret;

  for(int i=0;i<nn;i++){
    for(int k=0;k<nn;k++){
      if(k != i){
        for(int j=0;j<pp;j++){
          tmpvec(j)=0.0;
        }
        tmp=0.0;

        for(int j=0;j<pp;j++){
          tmpvec(j)+=des(i*pp+j)-des(k*pp+j);
          tmp +=pow(tmpvec(j),2);
        }
        tmpvec=tmpvec / pow(tmp,(2.0-q)/2.0);

        for(int j=0;j<pp;j++){
          ret2(i*pp+j)+=tmpvec(j);
        }
      }
    }
  }
  ret2=(2.0 * q)/((double) nn * nn)* ret2;

  return(ret-ret2);
}

// // [[Rcpp::export]]

// // [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
NumericMatrix sp_cpp(int des_num,int dim_num,NumericMatrix& ini,
                            NumericVector& distind,List distparam,
                            NumericMatrix& distsamp,bool thin,NumericMatrix& bd,
                            int point_num,int it_max,int it_min,double tol,
                            int num_proc,double n0,NumericVector& wts,bool rnd_flg){

  int it_num=0;
  bool cont=true;
  bool cont2=true;
  arma::vec curconst(des_num);
  curconst.fill(0.0);
  arma::vec runconst(des_num);
  runconst.fill(0.0);
  arma::vec runconst_up(des_num);
  runconst_up.fill(0.0);

  std::vector<double> prevdes(des_num*dim_num);
  std::default_random_engine generator;
  generator.seed(std::random_device{}());
  int distint=0;

  std::vector<double> des(des_num*dim_num);
  std::vector<double> des_up(des_num*dim_num);
  for(int i=0;i<des_num;i++){
    for(int j=0;j<dim_num;j++){
      des[j+i*dim_num]=ini(i,j);
    }
  }

  double nug=0.0;
  Rcout << "Optimizing ... " << endl;
  while(cont){

    double prop=(double) it_num/(double) it_max;
    if(it_num>0){printBar(prop);}

    curconst.fill(0.0);
    bool nanflg=false;

    for(int i=0;i<des_num;i++){
      for(int j=0;j<dim_num;j++){
        prevdes[j+i*dim_num]=des[j+i*dim_num];
      }
    }

    arma::mat rnd(point_num,dim_num);
    arma::vec rnd_wts(point_num);
    if(thin){
      std::uniform_int_distribution<int> uddist(0,distsamp.nrow()-1);
      for(int i=0;i<point_num;i++){
        int ss;
        if(rnd_flg){
          ss=uddist(generator);
        }else{
          ss=i;
        }
        for(int j=0;j<dim_num;j++){
          rnd(i,j)=distsamp(ss,j);
        }
        rnd_wts(i)=wts(ss);
      }
    }else{
      for(int n=0;n<dim_num;n++){

        for(int i=0;i<point_num;i++){
          rnd_wts(i)=wts(i);
        }

        distint=distind(n);
        switch(distint){
        case 1:
        {
          SEXP unill=distparam[n];
          NumericVector uniyy(unill);
          std::uniform_real_distribution<double> uni_dist(uniyy(0),uniyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=uni_dist(generator);
          }
          break;
        }
        case 2:
        {
          SEXP normll=distparam[n];
          NumericVector normyy(normll);
          std::normal_distribution<double> norm_dist(normyy(0),normyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=norm_dist(generator);
          }
          break;
        }
        case 3:
        {
          SEXP expll=distparam[n];
          NumericVector expyy(expll);
          std::exponential_distribution<double> exp_dist(expyy(0));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=exp_dist(generator);
          }
          break;
        }
        case 4:
        {
          SEXP gamll=distparam[n];
          NumericVector gamyy(gamll);
          std::gamma_distribution<double> gam_dist(gamyy(0),gamyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=gam_dist(generator);
          }
          break;
        }
        case 5:
        {
          SEXP lnll=distparam[n];
          NumericVector lnyy(lnll);
          std::lognormal_distribution<double> ln_dist(lnyy(0),lnyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=ln_dist(generator)/(exp(pow(lnyy(1),2.0)-1.0)* exp(2*lnyy(0)+pow(lnyy(1),2.0)));
          }
          break;
        }
        case 6:
        {
          SEXP tll=distparam[n];
          NumericVector tyy(tll);
          std::student_t_distribution<double> t_dist(tyy(0));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=t_dist(generator);
          }
          break;
        }
        case 7:
        {
          SEXP wbll=distparam[n];
          NumericVector wbyy(wbll);
          std::weibull_distribution<double> wb_dist(wbyy(0),wbyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=wb_dist(generator);
          }
          break;
        }
        case 8:
        {
          SEXP cauchll=distparam[n];
          NumericVector cauchyy(cauchll);
          std::cauchy_distribution<double> cauch_dist(cauchyy(0),cauchyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=cauch_dist(generator);
          }
          break;
        }
        case 9:
        {
          SEXP betall=distparam[n];
          NumericVector betayy(betall);
          sftrabbit::beta_distribution<> beta_dist(betayy(0),betayy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=beta_dist(generator);
          }
          break;
        }
        }
      }
    }

    for(int m=0;m<des_num;m++){

      arma::vec xprime(dim_num);
      xprime.fill(0.0);
      arma::vec tmpvec(dim_num);
      tmpvec.fill(0.0);

      for(int o=0;o<des_num;o++){

        if(o != m){
          double tmptol=0.0;
          for(int n=0;n<dim_num;n++){
            tmpvec(n)=prevdes[n+m*dim_num]- prevdes[n+o*dim_num];
            tmptol +=pow(tmpvec(n),2.0);
          }
          tmptol=sqrt(tmptol);

          for(int n=0;n<dim_num;n++){
            xprime(n)+=tmpvec(n)/(tmptol+nug*DBL_MIN);
          }
        }
      }

      for(int n=0;n<dim_num;n++){
        xprime(n)=xprime(n)*((double) point_num/(double) des_num);
      }

      for(int o=0;o<point_num;o++){
        double tmptol=0.0;
        for(int n=0;n<dim_num;n++){
          tmptol +=pow(rnd(o,n)- prevdes[n+m*dim_num],2.0);
        }
        tmptol=sqrt(tmptol);
        curconst(m)+=rnd_wts(o)/(tmptol+(nug*DBL_MIN));

        for(int n=0;n<dim_num;n++){
          xprime(n)+=rnd_wts(o)*rnd(o,n)/(tmptol+(nug*DBL_MIN));
        }
      }

      double denom=(1.0-(n0/(it_num+n0)))* runconst(m)+(n0/(it_num+n0))* curconst(m);
      for(int n=0;n<dim_num;n++){
        xprime(n)=((1.0-(n0/(it_num+n0)))* runconst(m)* prevdes[n+m*dim_num]+(n0/(it_num+n0))* xprime(n))/ denom;
      }

      for(int n=0;n<dim_num;n++){
        xprime(n)=min(max(xprime(n),bd(n,0)),bd(n,1));
      }

      for(int n=0;n<dim_num;n++){
        des_up[n+m*dim_num]=xprime(n);
        if(isnan(xprime(n))){
          nanflg=true;
        }
      }
      runconst_up(m)=(1-(n0/(it_num+n0)))*runconst(m)+(n0/(it_num+n0))*curconst(m);

    }

    if(nanflg){
      nug +=1.0;
      runconst.fill(0.0);
      Rcout << endl << "Numerical instablities encountered... resetting optimization" << endl;
    }else{
      des=des_up;
      runconst=runconst_up;
    }

    it_num++;

    double maxdiff=0.0;
    double rundiff=0.0;
    for(int n=0;n<des_num;n++){
      rundiff=0.0;
      for(int o=0;o<dim_num;o++){
        rundiff +=std::pow(des[o+n*dim_num]-prevdes[o+n*dim_num],2.0);
      }
      maxdiff=max(maxdiff,rundiff);
    }
    if((maxdiff < tol)&&(it_num >= it_min)&&(!nanflg)){
      cont=false;
      Rcout << endl << "Tolerance level reached: done!" << endl;
    }

    if((it_num >= it_max)&&(!nanflg)){
      cont=false;
      Rcout << endl << "Maximum iterations reached: done!" << endl;
    }

  }

  NumericMatrix retdes(des_num,dim_num);
  for(int j=0;j<dim_num;j++){
      for(int i=0;i<des_num;i++){
        retdes(i,j)=des[j+i*dim_num];
      }
  }

  return(retdes);

}

// // [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
NumericMatrix sp_seq_cpp(NumericMatrix& cur,int nseq,NumericMatrix& ini,
                        NumericVector& distind,List distparam,
                        NumericMatrix& distsamp,bool thin,NumericMatrix& bd,
                        int point_num,int it_max,int it_min,double tol,int num_proc){

  int dim_num=cur.ncol();
  int nini=cur.nrow();

  std::default_random_engine generator;
  int distint=0;

  std::vector<double> des((nini+nseq)*dim_num);
  std::vector<double> des_up((nini+nseq)*dim_num);
  for(int i=0;i<nini;i++){
    for(int j=0;j<dim_num;j++){
      des[j+i*dim_num]=cur(i,j);
    }
  }
  for(int i=nini;i<(nini+nseq);i++){
    for(int j=0;j<dim_num;j++){
      des[j+i*dim_num]=ini(i-nini,j);
    }
  }
  des_up=des;

  bool cont=true;
  int it_cur=0;
  Rcout << "Optimizing ... " << endl;

  while(cont){

    double prop=(double) it_cur/(double) it_max;
    if(it_cur>0){printBar(prop);}

    rst:
      arma::mat rnd(point_num,dim_num);
    if(thin){
      std::uniform_int_distribution<int> uddist(0,distsamp.nrow()-1);
      for(int i=0;i<point_num;i++){
        int ss=uddist(generator);
        for(int j=0;j<dim_num;j++){
          rnd(i,j)=distsamp(ss,j);
        }
      }
    }else{
      for(int n=0;n<dim_num;n++){
        distint=distind(n);
        switch(distint){
        case 1:
        {
          SEXP unill=distparam[n];
          NumericVector uniyy(unill);
          std::uniform_real_distribution<double> uni_dist(uniyy(0),uniyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=uni_dist(generator);
          }
          break;
        }
        case 2:
        {
          SEXP normll=distparam[n];
          NumericVector normyy(normll);
          std::normal_distribution<double> norm_dist(normyy(0),normyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=norm_dist(generator);
          }
          break;
        }
        case 3:
        {
          SEXP expll=distparam[n];
          NumericVector expyy(expll);
          std::exponential_distribution<double> exp_dist(expyy(0));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=exp_dist(generator);
          }
          break;
        }
        case 4:
        {
          SEXP gamll=distparam[n];
          NumericVector gamyy(gamll);
          std::gamma_distribution<double> gam_dist(gamyy(0),gamyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=gam_dist(generator);
          }
          break;
        }
        case 5:
        {
          SEXP lnll=distparam[n];
          NumericVector lnyy(lnll);
          std::lognormal_distribution<double> ln_dist(lnyy(0),lnyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=ln_dist(generator)/(exp(pow(lnyy(1),2.0)-1.0)* exp(2*lnyy(0)+pow(lnyy(1),2.0)));
          }
          break;
        }
        case 6:
        {
          SEXP tll=distparam[n];
          NumericVector tyy(tll);
          std::student_t_distribution<double> t_dist(tyy(0));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=t_dist(generator);
          }
          break;
        }
        case 7:
        {
          SEXP wbll=distparam[n];
          NumericVector wbyy(wbll);
          std::weibull_distribution<double> wb_dist(wbyy(0),wbyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=wb_dist(generator);
          }
          break;
        }
        case 8:
        {
          SEXP cauchll=distparam[n];
          NumericVector cauchyy(cauchll);
          std::cauchy_distribution<double> cauch_dist(cauchyy(0),cauchyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=cauch_dist(generator);
          }
          break;
        }
        case 9:
        {
          SEXP betall=distparam[n];
          NumericVector betayy(betall);
          sftrabbit::beta_distribution<> beta_dist(betayy(0),betayy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=beta_dist(generator);
          }
          break;
        }
        }
      }
    }

    for(int m=0;m<nseq;m++){

      arma::vec xprime(dim_num);
      arma::vec prevvec(dim_num);
      for(int n=0;n<dim_num;n++){
        prevvec(n)=des[n+(nini+m)*dim_num];
      }
      xprime.fill(0.0);
      arma::vec tmpvec(dim_num);
      tmpvec.fill(0.0);

      double tmptol;
      for(int o=0;o<(nini+nseq);o++){
        tmptol=0.0;
        for(int n=0;n<dim_num;n++){
          tmpvec(n)=prevvec(n)- des[n+o*dim_num];
          tmptol +=pow(tmpvec(n),2.0);
        }
        tmptol=sqrt(tmptol);

        if(tmptol>0.0){
          for(int n=0;n<dim_num;n++){
            xprime(n)+=tmpvec(n)/tmptol;
          }
        }
      }

      for(int n=0;n<dim_num;n++){
        xprime(n)=xprime(n)*((double) point_num/(double)(nini+nseq));
      }

      double tmpconst=0.0;
      for(int o=0;o<point_num;o++){
        double tmptol=0.0;
        for(int n=0;n<dim_num;n++){
          tmpvec(n)=pow(rnd(o,n)- prevvec(n),2);
          tmptol +=tmpvec(n);
        }
        tmptol=sqrt(tmptol);
        tmpconst +=1.0/tmptol;

        for(int n=0;n<dim_num;n++){
          xprime(n)+=rnd(o,n)/tmptol;
        }
      }

      for(int n=0;n<dim_num;n++){
        xprime(n)=xprime(n)/tmpconst;
      }

      for(int n=0;n<dim_num;n++){
        xprime(n)=min(max(xprime(n),bd(n,0)),bd(n,1));
      }

      for(int n=0;n<dim_num;n++){
        des_up[n+(nini+m)*dim_num]=xprime(n);
      }

    }

    des=des_up;

    it_cur++;
    if(it_cur >= it_max){
      cont=false;
    }
  }

  NumericMatrix retdes(nini+nseq,dim_num);
  for(int j=0;j<dim_num;j++){
    for(int i=0;i<(nini+nseq);i++){
      retdes(i,j)=des[j+i*dim_num];
    }
  }

  return(retdes);

}

// [[Rcpp::export]]
double gamma_eval(arma::vec& dd,arma::vec& theta){
  double runtot=0.0;
  int pp=dd.size();
  for(int i=0;i<pp;i++){
    runtot +=theta[i]*pow(dd[i],2.0);
  }
  runtot=exp(-runtot);
  return(runtot);
}

arma::vec asvec(arma::rowvec xx){
  arma::vec ret(xx.size());
  for(int i=0;i<xx.size();i++){
    ret(i)=xx(i);
  }
  return(ret);
}

// [[Rcpp::export]]
arma::vec omega(arma::vec& theta_vec,arma::vec& gamma_vec,int max_ord){
  arma::vec theta_tr(theta_vec.size()-1);
  int ct;
  int len=std::min((int) theta_vec.size()-1,max_ord);
  arma::vec rec_vec(len+1);
  arma::vec ret(theta_vec.size());

  for(int i=0;i<theta_vec.size();i++){
    ct=0;
    for(int j=0;j<theta_vec.size();j++){
      if(j == i){
      }else{
        theta_tr(ct)=theta_vec(j);
        ct++;
      }
    }
    for(int j=0;j<(len+1);j++){
      if(j==0){
        rec_vec(j)=1.0;
      }
      else{
        rec_vec(j)=0.0;
      }
    }
    for(int j=0;j<(theta_vec.size()-1);j++){
      rec_vec(0)=rec_vec(0);
      rec_vec(arma::span(1,rec_vec.size()-1))=theta_tr(j)*rec_vec(arma::span(0,rec_vec.size()-2))+ rec_vec(arma::span(1,rec_vec.size()-1));
    }
    ret(i)=theta_vec(i)* dot(rec_vec(arma::span(0,max_ord-1)),gamma_vec(arma::span(0,max_ord-1)));

  }
  return(ret);
}

// [[Rcpp::export]]
arma::mat opt_hess(arma::vec zz,arma::vec omega_vec){
  int pp=zz.size();
  arma::mat retmtx(pp,pp);

  arma::vec zz_tilde=zz;
  double nm=0.0;
  for(int i=0;i<pp;i++){
    zz_tilde(i)=zz(i)*omega_vec(i);
    nm +=pow(zz(i)*omega_vec(i),2.0);
  }

  for(int i=0;i<pp;i++){
    for(int j=0;j<pp;j++){
      retmtx(i,j)=4.0/exp(1.0)/nm*zz_tilde(i)*zz_tilde(j);
    }
  }
  return(retmtx);
}

// [[Rcpp::export]]
arma::vec psp_mi(arma::vec& xx,arma::mat& omega_mat,arma::mat& samp_mat,
                std::vector<double>& des_mat,int des_num,int ii){

  int nn=des_num;
  int pp=samp_mat.n_cols;
  int NN=samp_mat.n_rows;
  int RR=omega_mat.n_rows;
  arma::vec pt_vec(pp);
  arma::vec omega_vec(pp);
  double gam=0.0;
  arma::vec runvec(pp);
  runvec.fill(0.0);
  arma::mat runmat(pp,pp);
  runmat.fill(0.0);
  arma::vec tmpvec(pp);
  arma::mat tmpmat(pp,pp);

  for(int r=0;r<RR;r++){
    omega_vec=omega_mat.col(r);
    for(int j=0;j<nn;j++){
      if(!(ii == j)){
        for(int l=0;l<pp;l++){
          pt_vec(l)=des_mat[l+j*pp];
        }
        pt_vec=xx - pt_vec;
        gam=gamma_eval(pt_vec,omega_vec);
        for(int l=0;l<pp;l++){
          runvec[l]+=(2.0)/(nn*RR)*gam*omega_vec(l)*pt_vec(l);
        }
        tmpmat=opt_hess(pt_vec,omega_vec);
        tmpvec=4.0/(nn*RR)*tmpmat*xx;
        runvec +=tmpvec;
      }
    }
  }
  for(int r=0;r<RR;r++){
    omega_vec=omega_mat.col(r);
    for(int m=0;m<NN;m++){
      pt_vec=asvec(samp_mat.row(m));
      pt_vec=xx - pt_vec;
      gam=gamma_eval(pt_vec,omega_vec);
      pt_vec=asvec(samp_mat.row(m));
      for(int l=0;l<pp;l++){
        runvec(l)+=(2.0)/(NN*RR)*gam*omega_vec(l)*pt_vec(l);
      }
    }
  }

  for(int r=0;r<RR;r++){
    omega_vec=omega_mat.col(r);
    for(int j=0;j<nn;j++){
      if(!(ii == j)){
        for(int l=0;l<pp;l++){
          pt_vec(l)=des_mat[l+j*pp];
        }
        pt_vec=xx - pt_vec;
        tmpmat=opt_hess(pt_vec,omega_vec);
        runmat +=4.0/(nn*RR)*tmpmat;
      }
    }
  }

  for(int r=0;r<RR;r++){
    omega_vec=omega_mat.col(r);
    for(int m=0;m<NN;m++){
      pt_vec=asvec(samp_mat.row(m));
      pt_vec=xx - pt_vec;
      gam=gamma_eval(pt_vec,omega_vec);
      for(int l=0;l<pp;l++){
        runmat(l,l)+=(2.0)/(NN*RR)*gam*omega_vec(l);
      }
    }
  }

  tmpmat=inv(runmat);
  tmpvec=tmpmat*runvec;
  return(tmpvec);
}

// [[Rcpp::export]]
NumericMatrix psp_cpp(NumericMatrix& Rcpp_inides,NumericVector& distind,List distparam,
                      NumericVector& gam_param,arma::vec& gamma_vec,int max_ord,
                      NumericMatrix& distsamp,bool thinind,
                      NumericMatrix& gamsamp,bool gamind,NumericMatrix& bd,
                      int point_num,int gam_point_num,int it_max,double tol,int num_proc){

  double gam_shape=gam_param[0];
  double gam_rate=gam_param[1];

  int des_num=Rcpp_inides.nrow();
  int dim_num=Rcpp_inides.ncol();

  int it_num=0;
  bool cont=true;
  bool cont2=true;

  std::vector<double> prevdes(des_num*dim_num);
  std::vector<double> tmpvec(dim_num);
  std::default_random_engine generator;
  int distint=0;

  std::vector<double> des(des_num*dim_num);
  for(int i=0;i<des_num;i++){
    for(int j=0;j<dim_num;j++){
      des[j+i*dim_num]=Rcpp_inides(i,j);
    }
  }

  arma::mat rnd(point_num,dim_num);
  arma::mat theta_mat(gam_point_num,dim_num);
  while(cont){

    arma::vec omega_vec(dim_num);
    arma::mat omega_mat(dim_num,gam_point_num);
    arma::vec theta_vec(dim_num);

    Rcout << "PSP: Iteration " << it_num << "/" << it_max << endl;

    for(int i=0;i<des_num;i++){
      for(int j=0;j<dim_num;j++){
        prevdes[j+i*dim_num]=des[j+i*dim_num];
      }
    }

    if(thinind){
      std::uniform_int_distribution<int> uddist(0,distsamp.nrow()-1);
      for(int i=0;i<point_num;i++){
        int ss=uddist(generator);
        for(int j=0;j<dim_num;j++){
          rnd(i,j)=distsamp(ss,j);
        }
      }
    }else{
      for(int n=0;n<dim_num;n++){
        distint=distind(n);
        switch(distint){
        case 1:
        {
          SEXP unill=distparam[n];
          NumericVector uniyy(unill);
          std::uniform_real_distribution<double> uni_dist(uniyy(0),uniyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=uni_dist(generator);
          }
          break;
        }
        case 2:
        {
          SEXP normll=distparam[n];
          NumericVector normyy(normll);
          std::normal_distribution<double> norm_dist(normyy(0),normyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=norm_dist(generator);
          }
          break;
        }
        case 3:
        {
          SEXP expll=distparam[n];
          NumericVector expyy(expll);
          std::exponential_distribution<double> exp_dist(expyy(0));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=exp_dist(generator);
          }
          break;
        }
        case 4:
        {
          SEXP gamll=distparam[n];
          NumericVector gamyy(gamll);
          std::gamma_distribution<double> gam_dist(gamyy(0),gamyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=gam_dist(generator);
          }
          break;
        }
        case 5:
        {
          SEXP lnll=distparam[n];
          NumericVector lnyy(lnll);
          std::lognormal_distribution<double> ln_dist(lnyy(0),lnyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=ln_dist(generator)/(exp(pow(lnyy(1),2.0)-1.0)* exp(2*lnyy(0)+pow(lnyy(1),2.0)));
          }
          break;
        }
        case 6:
        {
          SEXP tll=distparam[n];
          NumericVector tyy(tll);
          std::student_t_distribution<double> t_dist(tyy(0));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=t_dist(generator);
          }
          break;
        }
        case 7:
        {
          SEXP wbll=distparam[n];
          NumericVector wbyy(wbll);
          std::weibull_distribution<double> wb_dist(wbyy(0),wbyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=wb_dist(generator);
          }
          break;
        }
        case 8:
        {
          SEXP cauchll=distparam[n];
          NumericVector cauchyy(cauchll);
          std::cauchy_distribution<double> cauch_dist(cauchyy(0),cauchyy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=cauch_dist(generator);
          }
          break;
        }
        case 9:
        {
          SEXP betall=distparam[n];
          NumericVector betayy(betall);
          sftrabbit::beta_distribution<> beta_dist(betayy(0),betayy(1));
          for(int o=0;o<point_num;o++){
            rnd(o,n)=beta_dist(generator);
          }
          break;
        }
        }
      }
    }

    for(int m=0;m<des_num;m++){

      vec runMat(dim_num);
      vec runVec(dim_num);
      runMat.fill(0.0);
      runVec.fill(0.0);
      arma::vec xprime(dim_num);
      xprime.fill(0.0);
      arma::vec tmpvec(dim_num);
      tmpvec.fill(0.0);

      for(int n=0;n<dim_num;n++){
        xprime[n]=des[n+m*dim_num];
      }

      xprime=psp_mi(xprime,omega_mat,rnd,des,des_num,m);

      for(int n=0;n<dim_num;n++){
        xprime[n]=min(max(xprime[n],bd(n,0)),bd(n,1));
      }

      for(int n=0;n<dim_num;n++){
        des[n+m*dim_num]=xprime[n];
      }

    }

    it_num++;
    if(it_num >= it_max){
      cont=false;
    }

  }

  stop:
    NumericMatrix retdes(des_num,dim_num);
  for(int i=0;i<des_num;i++){
    for(int j=0;j<dim_num;j++){
      retdes(i,j)=des[j+i*dim_num];
    }
  }

  return(retdes);

}

// [[Rcpp::export]]
NumericMatrix psp_seq_cpp(NumericMatrix& cur,int nseq,NumericMatrix& ini,
                          NumericVector& distind,List distparam,
                          NumericVector& gam_param,arma::vec& gamma_vec,int max_ord,
                          NumericMatrix& distsamp,bool thin,
                          NumericMatrix& gamsamp,bool gamind,NumericMatrix& bd,
                          int point_num,int gam_point_num,int it_max,double tol,int num_proc){

  double gam_shape=gam_param[0];
  double gam_rate=gam_param[1];

  int dim_num=cur.ncol();
  int nini=cur.nrow();

  std::default_random_engine generator;
  int distint=0;

  std::vector<double> des((nini+nseq)*dim_num);
  for(int i=0;i<nini;i++){
    for(int j=0;j<dim_num;j++){
      des[j+i*dim_num]=cur(i,j);
    }
  }
  for(int i=nini;i<(nini+nseq);i++){
    for(int j=0;j<dim_num;j++){
      des[j+i*dim_num]=ini(i-nini,j);
    }
  }

  for(int m=0;m<nseq;m++){

    Rcout << "Sequential PSPs: point " <<(nini+m+1)<< endl;

    int it_cur=0;
    arma::vec xprime(dim_num);
    arma::vec prevvec(dim_num);
    for(int n=0;n<dim_num;n++){
      prevvec(n)=des[n+(nini+m)*dim_num];
    }
    bool cont=true;

    while(cont){
      arma::vec omega_vec(dim_num);
      arma::mat omega_mat(dim_num,gam_point_num);
      arma::vec theta_vec(dim_num);

      arma::mat rnd(point_num,dim_num);
      if(thin){
        std::uniform_int_distribution<int> uddist(0,distsamp.nrow()-1);
        for(int i=0;i<point_num;i++){
          int ss=uddist(generator);
          for(int j=0;j<dim_num;j++){
            rnd(i,j)=distsamp(ss,j);
          }
        }
      }else{
        for(int n=0;n<dim_num;n++){
          distint=distind(n);
          switch(distint){
          case 1:
          {
            SEXP unill=distparam[n];
            NumericVector uniyy(unill);
            std::uniform_real_distribution<double> uni_dist(uniyy(0),uniyy(1));
            for(int o=0;o<point_num;o++){
              rnd(o,n)=uni_dist(generator);
            }
            break;
          }
          case 2:
          {
            SEXP normll=distparam[n];
            NumericVector normyy(normll);
            std::normal_distribution<double> norm_dist(normyy(0),normyy(1));
            for(int o=0;o<point_num;o++){
              rnd(o,n)=norm_dist(generator);
            }
            break;
          }
          case 3:
          {
            SEXP expll=distparam[n];
            NumericVector expyy(expll);
            std::exponential_distribution<double> exp_dist(expyy(0));
            for(int o=0;o<point_num;o++){
              rnd(o,n)=exp_dist(generator);
            }
            break;
          }
          case 4:
          {
            SEXP gamll=distparam[n];
            NumericVector gamyy(gamll);
            std::gamma_distribution<double> gam_dist(gamyy(0),gamyy(1));
            for(int o=0;o<point_num;o++){
              rnd(o,n)=gam_dist(generator);
            }
            break;
          }
          case 5:
          {
            SEXP lnll=distparam[n];
            NumericVector lnyy(lnll);
            std::lognormal_distribution<double> ln_dist(lnyy(0),lnyy(1));
            for(int o=0;o<point_num;o++){
              rnd(o,n)=ln_dist(generator)/(exp(pow(lnyy(1),2.0)-1.0)* exp(2*lnyy(0)+pow(lnyy(1),2.0)));
            }
            break;
          }
          case 6:
          {
            SEXP tll=distparam[n];
            NumericVector tyy(tll);
            std::student_t_distribution<double> t_dist(tyy(0));
            for(int o=0;o<point_num;o++){
              rnd(o,n)=t_dist(generator);
            }
            break;
          }
          case 7:
          {
            SEXP wbll=distparam[n];
            NumericVector wbyy(wbll);
            std::weibull_distribution<double> wb_dist(wbyy(0),wbyy(1));
            for(int o=0;o<point_num;o++){
              rnd(o,n)=wb_dist(generator);
            }
            break;
          }
          case 8:
          {
            SEXP cauchll=distparam[n];
            NumericVector cauchyy(cauchll);
            std::cauchy_distribution<double> cauch_dist(cauchyy(0),cauchyy(1));
            for(int o=0;o<point_num;o++){
              rnd(o,n)=cauch_dist(generator);
            }
            break;
          }
          case 9:
          {
            SEXP betall=distparam[n];
            NumericVector betayy(betall);
            sftrabbit::beta_distribution<> beta_dist(betayy(0),betayy(1));
            for(int o=0;o<point_num;o++){
              rnd(o,n)=beta_dist(generator);
            }
            break;
          }
          }
        }
      }

      arma::mat theta_mat(gam_point_num,dim_num);
      if(gamind){
        std::uniform_int_distribution<int> uddist(0,gamsamp.nrow()-1);
        for(int i=0;i<gam_point_num;i++){
          int ss=uddist(generator);
          for(int j=0;j<dim_num;j++){
            theta_mat(i,j)=gamsamp(ss,j);
          }
        }
      }
      else{
        std::gamma_distribution<double> dist(gam_shape,1.0/gam_rate);
        for(int i=0;i<gam_point_num;i++){
          for(int n=0;n<dim_num;n++){
            theta_mat(i,n)=dist(generator);
          }
        }
      }

      for(int o=0;o<gam_point_num;o++){
        for(int n=0;n<dim_num;n++){
          theta_vec(n)=theta_mat(o,n);
        }
        omega_vec=omega(theta_vec,gamma_vec,max_ord);
        omega_mat.col(o)=omega_vec;
      }

      for(int n=0;n<dim_num;n++){
        xprime(n)=prevvec(n);
      }

      xprime=psp_mi(xprime,omega_mat,rnd,des,m,m);

      for(int n=0;n<dim_num;n++){
        xprime(n)=min(max(xprime(n),bd(n,0)),bd(n,1));
      }

      double rundiff=0.0;
      for(int n=0;n<dim_num;n++){
        rundiff +=abs(xprime(n)-prevvec(n));
      }
      if((rundiff < tol)||(it_cur>=it_max)){
        cont=false;
      }else{
        it_cur++;
      }

      for(int n=0;n<dim_num;n++){
        prevvec(n)=xprime(n);
      }
    }

    for(int n=0;n<dim_num;n++){
      des[n+(nini+m)*dim_num]=xprime(n);
    }
  }

  NumericMatrix retdes(nini+nseq,dim_num);
  for(int j=0;j<dim_num;j++){
    for(int i=0;i<(nini+nseq);i++){
      retdes(i,j)=des[j+i*dim_num];
    }
  }

  return(retdes);

}

// // [[Rcpp::export]]

// // [[Rcpp::export]]

// // [[Rcpp::export]]

// // [[Rcpp::export]]

// // // [[Rcpp::export]]