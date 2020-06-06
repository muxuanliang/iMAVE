

#include "utils.h"
#include <cmath>

double threshold(double num) 
{
  return num > 0 ? num : 0;
}

// computes cumulative sum of vector x
VectorXd cumsum(const VectorXd& x) {
  const int n(x.size());
  VectorXd cmsm(n);
  cmsm(0) = x(0);
  
  for (int i = 1; i < n; i++) {
    cmsm(i) = cmsm(i-1) + x(i);
  }
  return (cmsm);
}

// computes reverse cumulative sum of vector x
VectorXd cumsumrev(const VectorXd& x) {
  const int n(x.size());
  VectorXd cmsm(n);
  cmsm(0) = x(n-1);
  
  for (int i = 1; i < n; i++) {
    cmsm(i) = cmsm(i-1) + x(n-i-1);
  }
  std::reverse(cmsm.data(), cmsm.data() + cmsm.size());
  return (cmsm);
}

VectorXd Standardize(const VectorXd& xx) {
  VectorXd xx_trans;
  int size = xx.size();
  double mean = xx.sum()/size;
  xx_trans = xx.array() - mean;
  double sd = sqrt(xx_trans.array().square().sum()/size);
  xx_trans = xx_trans/sd;
  return xx_trans;
}

MatrixXd Stdm(const MatrixXd& xx) {
  MatrixXd xx_trans(xx.rows(),xx.cols());
  int size = xx.cols();
  for (int i = 0; i<= size-1; i++)
  {
    xx_trans.col(i) = Standardize(xx.col(i)); 
  }
  return xx_trans;
}

MatrixXd kernel_std(const MatrixXd& kernel) {
  MatrixXd new_kernel(kernel.rows(), kernel.cols());
  for (int i = 0; i <= kernel.cols()-1; i++)
  {
    new_kernel.col(i) = kernel.col(i) / kernel.col(i).sum();
  }
  return new_kernel;
}

VectorXd VV(const VectorXd& xx, const VectorXd& yy) {
  int dimx = xx.size();
  int dimy = yy.size();
  VectorXd vv(dimx * dimy);
  
  for (int i = 0; i <= dimx-1; i++)
  {
    vv.segment(i * dimy ,dimy) = xx(i) * yy;
  }
  
  return vv;
}


//computes X'WX where W is diagonal (input w as vector)
MatrixXd XtWX(const MatrixXd& xx, const MatrixXd& ww) {
  const int n(xx.cols());
  MatrixXd AtWA(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(xx.adjoint() * ww.array().sqrt().matrix().asDiagonal()));
  return (AtWA);
}

//computes X'WX where W is diagonal (input w as vector)
MatrixXd XWXt(const MatrixXd& xx, const MatrixXd& ww) {
  const int n(xx.rows());
  MatrixXd AWAt(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(xx * ww.array().sqrt().matrix().asDiagonal()));
  return (AWAt);
}

//computes X'X
MatrixXd XtX(const MatrixXd& xx) {
  const int n(xx.cols());
  MatrixXd AtA(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(xx.adjoint()));
  return (AtA);
}

//computes XX'
MatrixXd XXt(const MatrixXd& xx) {
  const int n(xx.rows());
  MatrixXd AAt(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(xx));
  return (AAt);
}

//computes X'X
SpMat XtX(const SpMat& xx) {
  const int n(xx.cols());
  SpMat AtA(SpMat(n, n).selfadjointView<Upper>().rankUpdate(xx.adjoint()));
  return (AtA);
}

//computes XX'
SpMat XXt(const SpMat& xx) {
  const int n(xx.rows());
  SpMat AAt(SpMat(n, n).selfadjointView<Upper>().rankUpdate(xx));
  return (AAt);
}

//computes X'WX where W is diagonal (input w as vector)
SpMat XtWX(const SpMat& xx, const MatrixXd& ww) {
  const int n(xx.cols());
  SpMat AtWA(SpMat(n, n).
    selfadjointView<Lower>().rankUpdate(xx.adjoint() * ww.array().sqrt().matrix().asDiagonal()));
  return (AtWA);
}

//computes X'WX where W is diagonal (input w as vector)
SpMat XWXt(const SpMat& xx, const MatrixXd& ww) {
  const int n(xx.rows());
  SpMat AWAt(SpMat(n, n).
    selfadjointView<Lower>().rankUpdate(xx * ww.array().sqrt().matrix().asDiagonal()));
  return (AWAt);
}


bool stopRule(const VectorXd& cur, const VectorXd& prev, const double& tolerance) {
  for (unsigned i = 0; i < cur.rows(); i++) {
    if ( (cur(i) != 0 && prev(i) == 0) || (cur(i) == 0 && prev(i) != 0) ) {
      return 0;
    }
    if (cur(i) != 0 && prev(i) != 0 && 
        std::abs( (cur(i) - prev(i)) / prev(i)) > tolerance) {
  	  return 0;
    }
  }
  return 1;
}

void createC(SpMatR &C, const SpMat& group, const int& M) {  

  int row_idx = 0;
  for (int k=0; k < group.outerSize(); ++k) 
  {
    for (SpMat::InnerIterator it(group, k); it; ++it)
    {
      C.insert(row_idx, it.row()) = 1.0;
      ++row_idx;
    }
  }
  
  C.makeCompressed();  
}
