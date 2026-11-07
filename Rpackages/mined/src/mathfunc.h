# ifndef MATHFUNCTIONS_H_
# define MATHFUNCTIONS_H_
# include <RcppEigen.h>
# include <unordered_set>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace RcppEigen;
using namespace Eigen;
using namespace std;

#define minTwo(a,b)(a<b? a:b)
#define MaxTwo(a,b)(a>b? a:b)

# define PI 3.1415926

bool isPrime(int n);

int reportMaxPrime(int n);

ComplexVector FFT(NumericVector& x);

ComplexVector IFFT(ComplexVector& x);

double calOmega(double x);

VectorXd calOmega(VectorXd & x);

NumericVector calOmega(NumericVector & x);

int POWMOD(int x,int a,int n);

int findPrimefactors(unordered_set<int> & s,int n);

VectorXd findFactorize(int n);
int generateOrp(int n);

int argMin(VectorXd & x);
int argMax(VectorXd & x);

VectorXd compareMin(VectorXd & x,VectorXd & y);

MatrixXd fastPdist(MatrixXd & x,MatrixXd & y);

VectorXd fastPdist(MatrixXd & x,VectorXd & y);

VectorXd quantileCPP(VectorXd & x,VectorXd & q);

MatrixXd varCPP(MatrixXd & x);
VectorXi sampleCPP(int dim);

MatrixXd subMatCols(MatrixXd & x,VectorXi & ID);

MatrixXd subMatRows(MatrixXd & x,VectorXi ID);

MatrixXd setDiagonal(MatrixXd & x,double y);

VectorXi orderCPP(VectorXd & x);
MatrixXd bindMatByRows(MatrixXd & x,MatrixXd & y);

MatrixXd removeMatRow(MatrixXd & x,int RowID);

MatrixXd bindMatByRows(MatrixXd & x,VectorXd & y);

template<typename Derived>
MatrixXd bindMatByRows(MatrixXd & x,Eigen::MatrixBase<Derived> & y);

MatrixXd removeMatCol(MatrixXd & x,int ColID);

MatrixXd removeMatRowCol(MatrixXd & x,int RowID,int ColID);

MatrixXd bindMatByCols(MatrixXd & x,MatrixXd & y);

MatrixXd bindMatByCols(MatrixXd & x,VectorXd & y);

template<typename Derived>
MatrixXd bindMatByCols(MatrixXd & x,Eigen::MatrixBase<Derived> & y);
VectorXd subVectElements(const VectorXd & x,const VectorXi EleID);

VectorXd bindVectorByRow(const VectorXd & x,const VectorXd & y);
# endif