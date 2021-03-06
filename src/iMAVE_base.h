#ifndef iMAVEBASE_H
#define iMAVEBASE_H


#include <RcppEigen.h>
#include <algorithm>
#include "utils.h"

template<typename VecTypeBeta, typename VecTypeAB, typename VecTypeKernel>
class iMAVEBase
{
protected:
  typedef double Double;
  typedef Eigen::Matrix<Double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
  
  const int nobs;
  const int nvars;
  const int dim;
  
  VecTypeBeta main_beta;
  VecTypeAB   main_AB;
  VecTypeKernel    main_kernel;
  VecTypeKernel    main_kernel_init;
  double      init_tol;
  double      tol;
  double      delta;
  
public:
  iMAVEBase(int n_, int p_, int d_, double tol_):
  nobs(n_), nvars(p_), dim(d_), init_tol(1e-3), tol(tol_), delta(999)
  {}
  
  virtual ~iMAVEBase(){}
  
  void update_AB(VecTypeKernel &kernel)
  {
    VecTypeAB newAB;
    next_AB(newAB, kernel);
    main_AB = newAB;
  }
  
  void update_beta(VecTypeKernel &kernel)
  {
    VecTypeBeta newbeta(nvars,dim);
    
    next_beta(newbeta, kernel);
    
    if (dim == 1)
    {
      delta = std::min((newbeta-main_beta).norm(),(newbeta+main_beta).norm());
    } else {
      delta = std::min((newbeta-main_beta).array().abs().maxCoeff(),(newbeta+main_beta).array().abs().maxCoeff());
    }
    
    main_beta = newbeta;
  }
  void update_kernel()
  {
    VecTypeKernel newkernel;
    next_kernel(newkernel);
    main_kernel = newkernel;
  }
  
  virtual void next_AB( VecTypeAB &res, const VecTypeKernel &kernel)=0;
  virtual void next_beta( VecTypeBeta &res, const VecTypeKernel& kernel)=0;
  virtual void next_kernel( VecTypeKernel& res)=0;
  virtual void init_kernel()=0;
  
  virtual void init_warm()
  {
    init_kernel();
    MatrixXd eye = Eigen::MatrixXd::Identity(nvars,nvars);
    if (dim == 1){
      main_beta = eye.col(0);
    } else {
      main_beta = eye.leftCols(dim);
    }
  }
  
  virtual int solve(int maxit){
    int i;
    for (i = 1, delta = 999; delta >= tol && i <= maxit; i++)
    {
      update_kernel();
      update_AB(main_kernel);
      update_beta(main_kernel);
    }
    return i;
  }
  
  
  virtual VecTypeBeta get_beta(){ return main_beta; }
  virtual VecTypeAB get_ab(){ return main_AB; }
  
};
#endif // iMAVEBASE_H
