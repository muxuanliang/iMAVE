#ifndef iMAVEZERO_H
#define iMAVEZERO_H

#include "iMAVE_base.h"
//#include "Linalg/BlasWrapper.h"
//#include "Spectra/SymEigsSolver.h"
#include <Eigen/Cholesky>
#include "utils.h"
#include <cmath>
#include <math.h>

using Rcpp::IntegerVector;

// MAVE and expansion inlcudes all data
// 
class iMAVEZERO: public iMAVEBase<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>
{
protected:
  typedef float Scalar;
  typedef double Double;
  typedef Eigen::Matrix<Double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
  typedef Eigen::Matrix<Double, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Map<const Matrix> MapMat;
  typedef Eigen::Map<const Vector> MapVec;
  typedef const Eigen::Ref<const Matrix> ConstGenericMatrix;
  typedef const Eigen::Ref<const Vector> ConstGenericVector;
  
  MapMat datX;                  // data matrix
  MapVec datY;
  MapVec datTr;
  MapVec datPi;
  
  int constraint;
  
  //using iMAVEBase <Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>:: nobs;
  
  
  void next_AB(VectorXd& res, const VectorXd& kernel)
  {
    VectorXd Y_trans(nobs);
    MatrixXd X_trans = Eigen::MatrixXd::Ones(nobs,dim+1) * 0.25;
    res.resize(2);
    
    Y_trans = datY.array() * (datTr.array() - 1/2);
    
    VectorXd link(nobs);
    link = datX * main_beta;
    
    X_trans.col(1) = link * 0.25;
    if (constraint == 0){
      res = weighted_least_sq(X_trans, Y_trans, kernel.array() * (datPi.array() * (datTr.array())+(1-datPi.array()) * (1-datTr.array())).array().inverse());
    } else {
      res = weighted_least_sqct(X_trans, Y_trans, kernel.array() * (datPi.array() * (datTr.array())+(1-datPi.array()) * (1-datTr.array())).array().inverse());
    }
  }
  
  void next_beta(Vector &res, const VectorXd& kernel)
  {
    VectorXd weight(nobs);
    VectorXd Y_trans(nobs);
    MatrixXd X_trans(nobs,nvars);
    VectorXd beta_resize(nvars);
    res.resize(nvars);
    
    for (int i = 0; i <= nobs-1; i++)
    {
      weight(i) = kernel(i) / (datPi(i) * datTr(i)+(1-datPi(i)) * (1-datTr(i)));
      Y_trans(i) = datY(i) * (datTr(i)-0.5) - 0.25 * main_AB(0);
      X_trans.row(i) = datX.row(i) * main_AB(1) * 0.25;
    }
    
    beta_resize = weighted_least_sq(X_trans, Y_trans, weight);
    
    double scale;
    scale = beta_resize.norm() * beta_resize(0) / fabs(beta_resize(0));
    main_AB(1) = main_AB(1) * scale;
    res = beta_resize/scale;
  }
  
  void next_kernel(VectorXd& res)
  {
    VectorXd link(nobs);
    link = datX * main_beta;
    res.resize(nobs);
    
    double h = pow(4/(dim+2),1/(dim+4)) * pow(nobs, -1/(dim+4));
    
    for (int i = 0; i <= nobs-1; i++)
    {
      res(i) = exp(-0.5 * pow(link(i),2) /(pow(h,2))) * pow(h,dim);
    }
    
    res = res/res.sum();
  }
  
  
public:
  iMAVEZERO(ConstGenericMatrix &datX_,
             ConstGenericVector &datY_,
             ConstGenericVector &datTr_,
             ConstGenericVector &datPi_,
             double tol_, int constraint_):
  iMAVEBase<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>(datX_.rows(), datX_.cols(), 1, tol_),
  datX(datX_.data(), datX_.rows(), datX_.cols()),
  datY(datY_.data(), datY_.size()),
  datTr(datTr_.data(), datTr_.size()),
  datPi(datPi_.data(), datPi_.size()),
  constraint(constraint_)
  {}
  
  virtual void init_kernel(){
    main_kernel_init.resize(nobs);
    
    double h = pow(4/(nvars+2),1/(nvars+4)) * pow(nobs, -1/(nvars+4));
    
    for (int i = 0; i <= nobs-1; i++)
    {
      main_kernel_init(i) = exp(-0.5*(datX.row(i).dot(datX.row(i)))/(pow(h,2))) * pow(h,nvars);
    }
    
    main_kernel_init = main_kernel_init/main_kernel_init.sum();
  }
  
  VectorXd weighted_least_sq(const MatrixXd & xx, const VectorXd & yy, const VectorXd & weight){
    MatrixXd XTWX = xx.adjoint() * weight.asDiagonal() * xx;
    VectorXd XTWY = xx.adjoint() * weight.asDiagonal() * yy;
    VectorXd AB = XTWX.llt().solve(XTWY);
    return AB;
  }
  
  VectorXd weighted_least_sqct(const MatrixXd & xx, const VectorXd & yy, const VectorXd & weight){
    MatrixXd XTWX = xx.adjoint() * weight.asDiagonal() * xx;
    VectorXd XTWY = xx.adjoint() * weight.asDiagonal() * yy;
    VectorXd AB = XTWX.llt().solve(XTWY);
    VectorXd new_AB = VectorXd::Zero(dim + 1);
    VectorXd difVec = VectorXd::Zero(dim + 1);
    VectorXd Grad = VectorXd::Zero(dim + 1);
    double step = 0.1;
    double dif = 1;
    int maxit = 1e5;
    for (int iter = 0; iter <= maxit & dif > 1e-5; iter++)
    {
      Grad = -xx.adjoint() * weight.asDiagonal() * (yy - xx * AB);
      new_AB = AB - step * Grad;
      new_AB(1) = new_AB(1) * (new_AB(1)>0);
      difVec = new_AB-AB;
      dif = difVec.array().abs().maxCoeff();
      AB = new_AB;
    }
    
    return AB.transpose();
  }
};

#endif // iMAVEZERO_H
