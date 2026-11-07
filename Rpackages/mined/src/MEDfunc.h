# ifndef MEDFUNC_H_
# define MEDFUNC_H_
# include "mathfunc.h"
# include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace RcppEigen;
using namespace Eigen;
using namespace std;

double getLogGenDis(const VectorXd & x,const VectorXd & y,const double s)
{
  VectorXd temp=x - y;
  try
  {
    if(x.size()!= y.size())
    {
      throw std::length_error("vectors must have the same length.");
    }
    if(s > 1e-6)
    {
      temp=temp.array().abs().pow(s);
      return log(temp.sum())/ s;
    }
    else
    {
      temp=temp.array().abs().log()/ x.size();
      return temp.sum();
    }
  }catch(std::exception &ex)
  {
    forward_exception_to_r(ex);
    return -1.0;
  }
}

VectorXd getLogGenDis(const MatrixXd & x,const VectorXd & y,const double s)
{
  VectorXd res(x.rows());
  for(int i=0;i < res.size();i++)
  {
    res(i)=-1.0;
  }
  int n(y.size());
  try
  {
    if(x.cols()!= n)
    {
      throw std::length_error("The columns of x must be same as the size of vector y");
    }
    VectorXd temp(x.cols());
    for(int i=0;i < res.size();i++)
    {
      temp=x.row(i);
      res(i)=getLogGenDis(temp,y,s);
    }
    return res;
  }catch(std::exception &ex)
  {
    forward_exception_to_r(ex);
    return res;
  }
}

MatrixXd translateMatrix(MatrixXd & SigmaMatrix,MatrixXd & x)
{
  MatrixXd res(x.rows(),x.cols());
  res=x * SigmaMatrix;
  return res;
}

VectorXd translateVector(MatrixXd & SigmaMatrix,VectorXd & x)
{
  return x.transpose()* SigmaMatrix;
}

MatrixXd sqrtInvMatrix(MatrixXd & SigmaMatrix)
{
  SelfAdjointEigenSolver<MatrixXd> res(SigmaMatrix);
  VectorXd eigenvalues=res.eigenvalues();
  eigenvalues=1.0 / eigenvalues.array().sqrt();
  return res.eigenvectors()* eigenvalues.asDiagonal()* res.eigenvectors().transpose();
}

List chooseMED(MatrixXd & candidates,
                VectorXd & candlf,
                int n,
                MatrixXd & Sigmak,
                double gamma,
                double s
)
{
  int N(candidates.rows());
  int MaxIndex(0);
  int dim=candidates.cols();
  MatrixXd MEDPoints(n,candidates.cols());
  VectorXd lfMED(n);
  MatrixXd SqrtSigmaInv=sqrtInvMatrix(Sigmak);
  VectorXd tempPoint(candidates.cols());
  VectorXd tempCrit(N);
  VectorXd rowMinCrit=MatrixXd::Constant(N,1,1e16);
  double templf(0.0);

  MaxIndex=argMax(candlf);
  tempPoint=candidates.row(MaxIndex);
  templf=candlf(MaxIndex);
  MEDPoints.row(0)=tempPoint;
  lfMED(0)=templf;

  MatrixXd transCand=translateMatrix(SqrtSigmaInv,candidates);
  VectorXd transtempPoint=translateVector(SqrtSigmaInv,tempPoint);
  MatrixXd CritCandMEDMatrix=MatrixXd::Constant(N,n,1e16);
  for(int i=1;i < n;i++)
  {
    tempCrit=0.5 * gamma *(templf + candlf.array())+ dim * getLogGenDis(transCand,transtempPoint,s).array();
    rowMinCrit=compareMin(rowMinCrit,tempCrit);
    CritCandMEDMatrix.col(i)=tempCrit;
    MaxIndex=argMax(rowMinCrit);
    tempPoint=candidates.row(MaxIndex);
    templf=candlf(MaxIndex);
    MEDPoints.row(i)=tempPoint;
    lfMED(i)=templf;
    transtempPoint=translateVector(SqrtSigmaInv,tempPoint);
  }
  return List::create(Named("points")=MEDPoints,Named("logf")=lfMED);
}
# endif