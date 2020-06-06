#define EIGEN_DONT_PARALLELIZE

#include "iMAVE_base.h"
#include "iMAVE_zero.h"
#include "iMAVE_all.h"

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//
using Eigen::MatrixXf;
using Eigen::MatrixXd;
using Eigen::VectorXf;
using Eigen::VectorXd;
using Eigen::ArrayXf;
using Eigen::ArrayXd;
using Eigen::ArrayXXf;
using Eigen::Map;

using Rcpp::wrap;
using Rcpp::as;
using Rcpp::List;
using Rcpp::Named;
using Rcpp::CharacterVector;

typedef Map<VectorXd> MapVecd;
typedef Map<Eigen::MatrixXd> MapMatd;


RcppExport SEXP iMAVE_cpp(SEXP x_,
                           SEXP y_,
                           SEXP tr_,
                           SEXP d_,
                           SEXP pi_,
                           SEXP init_,
                           SEXP method_,
                           SEXP opts_
                           ) 
{
  BEGIN_RCPP
  
  Rcpp::NumericMatrix xx(x_);
  Rcpp::NumericVector yy(y_);
  Rcpp::NumericVector tr(tr_);
  Rcpp::NumericVector pi(pi_);
  
  const int n = xx.rows();
  const int p = xx.cols();
  
  MatrixXd datX(n, p);
  VectorXd datY(n);
  VectorXd datTr(n);
  VectorXd datPi(n);
  
  // Copy data and convert type from double to float
  std::copy(xx.begin(), xx.end(), datX.data());
  std::copy(yy.begin(), yy.end(), datY.data());
  std::copy(tr.begin(), tr.end(), datTr.data());
  std::copy(pi.begin(), pi.end(), datPi.data());
  
  List opts(opts_);
  const int maxit          = as<int>(opts["maxit"]);
  const double tol         = as<double>(opts["tol"]);
  const int d              = as<int>(d_);
  const int constraint     = as<int>(opts["constraint"]);
  const bool normalizeWeight = as<bool>(opts["normalizeWeight"]);
  const 
  CharacterVector method(as<CharacterVector>(method_));
  CharacterVector initial(as<CharacterVector>(init_));
  
   if (method(0) == "zero")
   {
  iMAVEBase<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> *solver = NULL;
  solver = new iMAVEZERO(datX, datY, datTr, datPi, tol, constraint);
  
  solver->init_warm();
  int itertmp = solver->solve(maxit);
  MatrixXd beta = solver->get_beta();
  MatrixXd ab = solver->get_ab();
  
  return List::create(Named("dim")      = d,
                      Named("beta")     = beta, 
                      Named("ab")       = ab,
                      Named("niter")    = itertmp);
   }
  else
  {
    iMAVEBase<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> *solver = NULL;
    if (initial(0) == "direct")
    {
      solver = new iMAVEAll(datX, datY, datTr, d, 0, datPi, tol, constraint, normalizeWeight);
    } else if (initial(0) == "zero")
    {
      solver = new iMAVEAll(datX, datY, datTr, d, 1, datPi, tol, constraint, normalizeWeight);
    }
    
    
    solver->init_warm();
    int itertmp = solver->solve(maxit * d);
    MatrixXd beta = solver->get_beta();
    MatrixXd ab = solver->get_ab();
    
    if (itertmp >= maxit && d==1)
    {
      iMAVEBase<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> *solver_zero = NULL;
      solver_zero = new iMAVEZERO(datX, datY, datTr, datPi, tol, 0);
      
      solver_zero->init_warm();
      solver_zero->solve(maxit);
      beta = solver_zero->get_beta();
      ab = solver_zero->get_ab();
      itertmp = 0;
    }

    
    return List::create(Named("dim")      = d,
                        Named("beta")     = beta, 
                        Named("ab")       = ab,
                        Named("niter")     =itertmp);

   }

  
 
  END_RCPP
}



