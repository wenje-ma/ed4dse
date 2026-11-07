# ifndef LKPREDICTOR_H_
# define LKPREDICTOR_H_
# include "mathfunc.h"
// [[Rcpp::plugins(cpp11)]]

class LKPredictor
{
public:
  LKPredictor(){};
  LKPredictor(MatrixXd & DesignPoints,VectorXd & LfVector)
  {
    m_DesignPoints=DesignPoints;
    m_LfVector=LfVector;
    this->N=DesignPoints.rows();
    this->Dim=DesignPoints.cols();
  }
  void setData(MatrixXd & DesignPoints,VectorXd & LfVector)
  {
    m_DesignPoints=DesignPoints;
    m_LfVector=LfVector;
    this->N=DesignPoints.rows();
    this->Dim=DesignPoints.cols();
  }
  void fit()
  {
    MatrixXd res=fastPdist(this->m_DesignPoints,this->m_DesignPoints);
    MatrixXd res_temp(res);
    double quant=MaxTwo(1.0 - 5.0 * this->Dim / this->N,.5);
    res_temp=res + MatrixXd::Identity(N,N)* 10.0 * this->Dim;
    VectorXd tempVector=res_temp.rowwise().minCoeff();
    VectorXd quantvect(1);
    quantvect << quant;
    this->dbar=2.0 * this->Dim * quant * quantileCPP(tempVector,quantvect)(0);
    res=1.0 /(1.0 + res.array().square()/(this->dbar * this->dbar));
    res_temp=res + MatrixXd::Identity(N,N)* 1e-6;
    this->corMatrixInv=res_temp.inverse();
    this->coef=this->corMatrixInv * this->m_LfVector;
    this->dnom=this->corMatrixInv.rowwise().sum();
  }
  List predictor(VectorXd & x)
  {
    VectorXd distVect=fastPdist(this->m_DesignPoints,x);
    double mind=distVect.minCoeff();
    distVect=1.0 /(1.0 + distVect.array().square()/(this->dbar * this->dbar));
    double pred=(distVect.transpose()* this->coef)[0]/(distVect.transpose()* this->dnom)[0];
    return List::create(_["mind"]=mind,_["pred"]=pred,Named("coef")=this->coef,Named("dnom")=this->dnom,
                        Named("corMatrixInv")=this->corMatrixInv);
  }
public:
  MatrixXd m_DesignPoints;
  VectorXd m_LfVector;
private:
  MatrixXd corMatrixInv;
  VectorXd coef;
  VectorXd dnom;
  double dbar;
  int N;
  int Dim;
};

# endif