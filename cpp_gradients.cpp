#include <Rcpp.h>
using namespace Rcpp;

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;

// [[Rcpp::export]]
List cpp_gradient_ab( const Eigen::ArrayXXd & Xs,
                  const Eigen::ArrayXd & Ys,
                  double sigma,
                  const std::vector<double> & alpha,
                  const std::vector<double> & gamma,
                  const std::vector<double> & a,
                  const std::vector<double> & b
){
  unsigned int n = Ys.size();
  double one_over_n = 1./n;
  double inv_sigma = 1./sigma;
  
  Eigen::MatrixXd K(n,n);
  Eigen::MatrixXd A(n,n);
  Eigen::MatrixXd B(n,n);
  
  double temp0 = alpha[0] + gamma[0], temp1 = alpha[1] + gamma[1];
  
  for(unsigned int i = 0 ; i < n ; i++){
    double tempi = (Xs(i,0)*temp0 + Xs(i,1)*temp1);
    for(unsigned int j = 0 ; j < n ; j++){
      K(i,j) = tempi*Ys(j);
      A(i,j) = a[i];
      B(i,j) = b[i];
    }
  }
  
  Eigen::MatrixXd Pi = Eigen::exp((K-A-B.transpose()).array()*inv_sigma);
  Eigen::MatrixXd Pitilde = Pi;
  Pitilde.row(0) = 0*Pi.row(0);
    
  Eigen::MatrixXd baseline1 = (Pi.array() * (Xs.col(0).matrix()*Ys.matrix().transpose()).array()).matrix();
  Eigen::MatrixXd baseline2 = (Pi.array() * (Xs.col(1).matrix()*Ys.matrix().transpose()).array()).matrix();
  
  Eigen::MatrixXd E(n,2);
  Eigen::MatrixXd F(2,n);
  
  E << baseline1.rowwise().sum(), baseline2.rowwise().sum();
  E.row(0) = 0*E.row(0);
  F << baseline1.colwise().sum(), baseline2.colwise().sum();
  
  Eigen::VectorXd diag_vec(Eigen::VectorXd::Ones(n));
  Eigen::MatrixXd diag = diag_vec.asDiagonal()*one_over_n;
  
  Eigen::MatrixXd LHS(2*n,2*n);
  LHS << diag, Pitilde,
           Pi.transpose(), diag;
  
  Eigen::MatrixXd RHS(2*n,2);
  RHS << E,
         F.transpose();
  
  Eigen::MatrixXd DaDb(2*n,2);
  DaDb = LHS.colPivHouseholderQr().solve(RHS);
  
  Eigen::MatrixXd Da(n,2);
  Eigen::MatrixXd Db(n,2);
  
  Da << DaDb.topRows(n);
  Db << DaDb.bottomRows(n);
    
  return List::create(Named("DaDb")=DaDb, Named("Da")=Da, Named("Db")=Db);
}
