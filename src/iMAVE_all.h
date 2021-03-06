#ifndef iMAVEALL_H
#define iMAVEALL_H

#include "iMAVE_base.h"
#include "utils.h"
#include <cmath>
#include <math.h>
#include <Eigen/Cholesky>

//using Rcpp::IntegerVector;

// MAVE and expansion inlcudes all data
// 
class iMAVEAll: public iMAVEBase<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>
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
  
  int    method_init;
  int    constraint;
  bool   normalizeWeight;
 
  
  void next_AB(Matrix &res, const Matrix& kernel)
  {
    res.setZero(nobs, dim+1);
    VectorXd Y_trans = datY.array() * (datTr.array() - 0.5);
  
    MatrixXd Xj_trans = Eigen::MatrixXd::Ones(nobs, dim+1) * 0.25;
    MatrixXd link =  datX * main_beta;
    
    for (int j = 0; j <= nobs-1; j++){
      if (dim == 1)
      {
        Xj_trans.col(1) = (link - link(j) * Eigen::VectorXd::Ones(nobs)) * 0.25;
      } else {
        for (int i = 1; i <= dim ; i++)
        {
          Xj_trans.col(i) = (link.col(i-1) - link(j,i-1) * Eigen::VectorXd::Ones(nobs)) * 0.25;
        }
      }
      if (constraint == 0){
        res.row(j) = weighted_least_sq(Xj_trans, Y_trans, kernel.col(j).array() * (datPi.array() * (datTr.array())+(1-datPi.array()) * (1-datTr.array())).array().inverse());
      } else {
        res.row(j) = weighted_least_sqct(Xj_trans, Y_trans, kernel.col(j).array() * (datPi.array() * (datTr.array())+(1-datPi.array()) * (1-datTr.array())).array().inverse());
      }
    }
  }
  
  void next_beta(MatrixXd &res, const Matrix& kernel)
  {
    res.setZero(nvars, dim);
    VectorXd weight(nobs * nobs);
    VectorXd Y_trans(nobs * nobs);
    MatrixXd X_trans(nobs * nobs,dim * nvars);
    VectorXd beta_resize(dim * nvars);
    
    VectorXd temp(dim+1);
    for (int i = 0; i <= nobs-1; i++)
    {
      for (int j = 0; j <= nobs-1; j++)
      {
        if (!normalizeWeight)
        {
          weight(i * nobs + j) = kernel(i,j) / (datPi(i) * datTr(i)+(1-datPi(i)) * (1-datTr(i)));
        } else {
          weight(i * nobs + j) = (kernel(i,j) / kernel.col(j).array().sum()) / (datPi(i) * datTr(i)+(1-datPi(i)) * (1-datTr(i)));
        }
        
        Y_trans(i * nobs + j) = datY(i) * (datTr(i)-0.5) - 0.25 * main_AB(j,0);
        temp = main_AB.row(j);
        X_trans.row(i * nobs + j) = 0.25 * VV(temp.segment(1,dim), datX.row(i)-datX.row(j));
      }
    }
    beta_resize = weighted_least_sq(X_trans, Y_trans, weight);
    
    double scale;
    if (dim == 1)
    {
      scale = beta_resize.norm()* beta_resize(0) / fabs(beta_resize(0));
      main_AB.col(1) = main_AB.col(1) * scale;
    
      res = beta_resize/scale;
    } else {
      for (int i = 0; i <= dim-1; i++)
      {
        res.col(i) = beta_resize.segment(i * nvars, nvars);
        VectorXd temp = res.col(i);
        scale = res.col(i).norm() * temp(0) / fabs(temp(0));
        res.col(i) = res.col(i)/scale;
        main_AB.col(i+1) = main_AB.col(i+1) * scale;
      }
    }
  }
  
  void next_kernel(Matrix& res)
  {
    res.setZero(nobs,nobs);
    MatrixXd link = datX * main_beta;
    
    double h = pow(4/(dim+2.0),1/(dim+4.0)) * pow(nobs, -1/(dim+4.0));
    
    for (int i = 0; i <= nobs-1; i++)
    {
      for (int j = 0; j <= nobs-1; j++)
      {
        if (dim == 1)
        {
          res(i,j) = exp(-0.5 * pow(link(i)-link(j),2) /(pow(h,2))) * pow(h,dim);
        }
        else
        {
          res(i,j) = exp(-0.5 * (link.row(i)-link.row(j)).dot(link.row(i)-link.row(j)) /(pow(h,2))) * pow(h,dim);
        }
      }
    }
    res = res * res.rows() / res.sum();
  }
  
  void init_warm()
  {
    init_kernel();
    MatrixXd eye = Eigen::MatrixXd::Identity(nvars,nvars);
    if (dim == 1){
      main_beta = eye.col(0);
      main_beta = main_beta/main_beta(0);
    } else {
      main_beta = eye.leftCols(dim);
    }
    
    if (method_init == 1)
    {
      iMAVEBase<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> *solver_zero = NULL;
      solver_zero = new iMAVEZERO(datX, datY, datTr, datPi, tol, constraint);
      solver_zero->init_warm();
      solver_zero->solve(10);
      main_beta = solver_zero->get_beta();
    }
  }
  
  
  
public:
  iMAVEAll(ConstGenericMatrix &datX_,
            ConstGenericVector &datY_,
            ConstGenericVector &datTr_,
            int dim_,
            int method_init_,
            ConstGenericVector &datPi_,
            double tol_, int constraint_, bool normalizeWeight_):
  iMAVEBase<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>(datX_.rows(), datX_.cols(), dim_, tol_),
  datX(datX_.data(), datX_.rows(), datX_.cols()),
  datY(datY_.data(), datY_.size()),
  datTr(datTr_.data(), datTr_.size()),
  datPi(datPi_.data(), datPi_.size()),
  method_init(method_init_),
  constraint(constraint_),
  normalizeWeight(normalizeWeight_)
  {}
  
  virtual void init_kernel(){
    main_kernel_init.resize(nobs,nobs);
    double h = pow(4/(nvars+2),1/(nvars+4)) * pow(nobs, -1/(nvars+4));
    
    for (int i = 0; i <= nobs-1; i++)
    {
      for (int j = 0; j <= nobs-1; j++)
      {
        main_kernel_init(i,j) = exp(-0.5*(datX.row(i)-datX.row(j)).dot(datX.row(i)-datX.row(j))/(pow(h,2))) * pow(h,nvars);
      }
    }
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

#endif // iMAVEALL_H
