#include <Rcpp.h>
using namespace Rcpp;

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;

// [[Rcpp::export]]
Eigen::ArrayXXd generate_xs(unsigned int n, unsigned int m, long seed){
  Eigen::ArrayXXd to_return(n,m);
  
  for(unsigned int i = 0 ; i < m ; i++){
    for(unsigned int j = 0 ; j < n ; j++) {
      to_return(j,i) = R::rnorm(0,1);
    }
  }
  return to_return;
}

// [[Rcpp::export]]
List cpp_IPFP( const Eigen::ArrayXXd & Xs,
           const Eigen::ArrayXd & Ys,
           double sigma,
           const std::vector<double> & alpha,
           const std::vector<double> & gamma,
           unsigned long max_iter,
           double epsilon
){
  unsigned int n = Ys.size();
  double one_over_n = 1./n;
  double inv_sigma = 1./sigma;
  
  Eigen::VectorXd A;
  Eigen::VectorXd B;
  
  Eigen::MatrixXd K(n,n);
  
  double temp0 = alpha[0] + gamma[0], temp1 = alpha[1] + gamma[1];
  
  for(unsigned int i = 0 ; i < n ; i++){
    double tempi = (Xs(i,0)*temp0 + Xs(i,1)*temp1)*inv_sigma;
    for(unsigned int j = 0 ; j < n ; j++)
      K(i,j) = std::exp(tempi*Ys(j));
  }
  
  B = Eigen::VectorXd::Ones(n);
  
  unsigned long iter_ipfp = 0;
  bool stop_condition = false;
  while(!stop_condition){
    iter_ipfp +=1 ;
    A = (1./(n*(K*B).array()));
    B = 1./(n*(A.transpose()*K).array());
    
    Eigen::ArrayXd test = (A.array() * (K*B).array()) - one_over_n ;
    if(max_iter<iter_ipfp || test.abs().maxCoeff()<epsilon) {
      stop_condition = true;
    }
  }
  printf("iter: %lu \n", iter_ipfp);
  //std::cout << Ys << std::endl;
  
  double norm = A(0);
  double inv_norm = 1/norm;
  
  Eigen::VectorXd A_norm = inv_norm * A;
  Eigen::VectorXd B_norm = norm * B;
  
  Eigen::VectorXd a = -sigma * Eigen::log(A_norm.array());
  Eigen::VectorXd b = -sigma * Eigen::log(B_norm.array());
  
  return List::create(Named("a")=a, Named("b")=b);
}

Eigen::MatrixXd fill_matrix( const Eigen::VectorXd & vec,
                  unsigned int n
){
  Eigen::MatrixXd mat(n,n);
  
  for(unsigned int i = 0 ; i < n ; i++){
    for(unsigned int j = 0 ; j < n ; j++)
      mat(i,j) = vec[i];
  }
  return mat;
}

// [[Rcpp::export]]
List cpp_IPFP_lse( const Eigen::ArrayXXd & Xs,
               const Eigen::ArrayXd & Ys,
               double sigma,
               const std::vector<double> & alpha,
               const std::vector<double> & gamma,
               unsigned long max_iter,
               double epsilon
){
  unsigned int n = Ys.size();
  double one_over_n = 1./n;
  double inv_sigma = 1./sigma;
  
  Eigen::VectorXd base = one_over_n * Eigen::VectorXd::Ones(n);
  
  Eigen::VectorXd mu = -sigma * Eigen::log(base.array());
  Eigen::VectorXd nu = -sigma * Eigen::log(base.array());
  
  Eigen::MatrixXd K(n,n);
  
  double temp0 = alpha[0] + gamma[0], temp1 = alpha[1] + gamma[1];
  
  for(unsigned int i = 0 ; i < n ; i++){
    double tempi = (Xs(i,0)*temp0 + Xs(i,1)*temp1)*inv_sigma;
    for(unsigned int j = 0 ; j < n ; j++)
      K(i,j) = tempi*Ys(j);
  }
  
  Eigen::VectorXd v = Eigen::VectorXd::Zero(n);
  Eigen::MatrixXd V = fill_matrix(v,n);

  Eigen::VectorXd vstar;
  Eigen::MatrixXd Vstar;  
  
  Eigen::VectorXd u;
  Eigen::MatrixXd U;
  
  Eigen::VectorXd ustar;
  Eigen::MatrixXd Ustar;  
  
  Eigen::ArrayXd  test;
  
  unsigned long iter_ipfp = 0;
  bool stop_condition = false;
  
   while(!stop_condition){
     iter_ipfp +=1 ;
      
     vstar = (K-V.transpose()).rowwise().maxCoeff();
     Vstar = fill_matrix(vstar,n);
     
     u = mu.array() + vstar.array() + sigma * Eigen::log( (Eigen::exp( (K - V.transpose() - Vstar).array()*inv_sigma)).rowwise().sum() ).array();
     U = fill_matrix(u,n);
     
     ustar = (K-U).colwise().maxCoeff();
     Ustar = fill_matrix(ustar,n);
     
     v = nu.array() + ustar.array() + sigma * Eigen::log( (Eigen::exp( (K - U - Ustar.transpose()).array()*inv_sigma)).transpose().rowwise().sum() ).array();
     V = fill_matrix(v,n);
     
     test = (u.array() - mu.array() - sigma * Eigen::log( (Eigen::exp( (K - V.transpose()).array()*inv_sigma)).rowwise().sum() ).array()) ;
     if(max_iter<iter_ipfp || test.rowwise().sum().abs().maxCoeff()<epsilon) {
       stop_condition = true;
     }
    }

  printf("iter: %lu \n", iter_ipfp);
  //std::cout << Ys << std::endl;
  
  Eigen::VectorXd A = Eigen::exp(-inv_sigma*u.array());
  Eigen::VectorXd B = Eigen::exp(-inv_sigma*v.array());
  
  double norm = A(0);
  double inv_norm = 1/norm;
  
  Eigen::VectorXd A_norm = inv_norm * A;
  Eigen::VectorXd B_norm = norm * B;
  
  Eigen::VectorXd a = -sigma * Eigen::log(A_norm.array());
  Eigen::VectorXd b = -sigma * Eigen::log(B_norm.array());

  return List::create(Named("a")=a, Named("b")=b);
}

// [[Rcpp::export]]
List cpp_IPFP_AffMat( const Eigen::ArrayXXd & Xs,
                   const Eigen::ArrayXd & Ys,
                   double sigma,
                   const Eigen::MatrixXd & AffMat,
                   unsigned long max_iter,
                   double epsilon
){
  unsigned int n = Ys.size();
  double one_over_n = 1./n;
  double inv_sigma = 1./sigma;
  
  Eigen::VectorXd A;
  Eigen::VectorXd B;
  
  Eigen::MatrixXd K(n,n);
  Eigen::MatrixXd rowxi;
  Eigen::MatrixXd rowyj;
  Eigen::MatrixXd tempi;
  
  Eigen::ArrayXd test1;
  Eigen::ArrayXd test2;
  
  for(unsigned int i = 0 ; i < n ; i++){
    rowxi = Xs.row(i);
    tempi = inv_sigma*rowxi*AffMat;
    for(unsigned int j = 0 ; j < n ; j++){
      rowyj = Ys.row(j).transpose();
      K(i,j) = std::exp((tempi*rowyj)(0,0));
    }
  }
  
  B = Eigen::VectorXd::Ones(n);
  
  unsigned long iter_ipfp = 0;
  bool stop_condition = false;
  while(!stop_condition){
    iter_ipfp +=1 ;
    A = (1./(n*(K*B).array()));
    B = 1./(n*(A.transpose()*K).array());
    
    test1 = (A.array() * (K*B).array()) - one_over_n ;
    test2 = (B.array() * (A.transpose()*K).array()) - one_over_n ;
    if(max_iter<iter_ipfp || (test1.abs().maxCoeff()<epsilon && test2.abs().maxCoeff()<epsilon )) {
      stop_condition = true;
    }
  }
  
  printf("iter IPFP: %lu \n", iter_ipfp);
  //std::cout << Ys << std::endl;

  double norm = A(0);
  double inv_norm = 1/norm;

  Eigen::VectorXd A_norm = inv_norm * A;
  Eigen::VectorXd B_norm = norm * B;

  Eigen::VectorXd a = -sigma * Eigen::log(A_norm.array());
  Eigen::VectorXd b = -sigma * Eigen::log(B_norm.array());
  
  return List::create(Named("a")=a, Named("b")=b, Named("iterIpfp")=iter_ipfp, Named("test2")=test2,  Named("test1")=test1);
}

// [[Rcpp::export]]
List cpp_IPFP_AffMat_lse( const Eigen::ArrayXXd & Xs,
                          const Eigen::ArrayXd & Ys,
                          double sigma,
                          const Eigen::MatrixXd & AffMat,
                          unsigned long max_iter,
                          double epsilon
){
  unsigned int n = Ys.size();
  double one_over_n = 1./n;
  double inv_sigma = 1./sigma;
  
  Eigen::VectorXd base = one_over_n * Eigen::VectorXd::Ones(n);
  
  Eigen::VectorXd mu = -sigma * Eigen::log(base.array());
  Eigen::VectorXd nu = -sigma * Eigen::log(base.array());
  
  Eigen::MatrixXd K(n,n);
  Eigen::MatrixXd rowxi;
  Eigen::MatrixXd rowyj;
  Eigen::MatrixXd tempi;
  
  for(unsigned int i = 0 ; i < n ; i++){
    rowxi = Xs.row(i);
    tempi = inv_sigma*rowxi*AffMat;
    for(unsigned int j = 0 ; j < n ; j++){
      rowyj = Ys.row(j).transpose();
      K(i,j) = (tempi*rowyj)(0,0);
    }
  }
  
  Eigen::VectorXd v = Eigen::VectorXd::Zero(n);
  Eigen::MatrixXd V = fill_matrix(v,n);
  
  Eigen::VectorXd vstar;
  Eigen::MatrixXd Vstar;  
  
  Eigen::VectorXd u;
  Eigen::MatrixXd U;
  
  Eigen::VectorXd ustar;
  Eigen::MatrixXd Ustar;  
  
  Eigen::ArrayXd  test;
  
  unsigned long iter_ipfp = 0;
  bool stop_condition = false;
  
  while(!stop_condition){
    iter_ipfp +=1 ;
    
    vstar = (K-V.transpose()).rowwise().maxCoeff();
    Vstar = fill_matrix(vstar,n);
    
    u = mu.array() + vstar.array() + sigma * Eigen::log( (Eigen::exp( (K - V.transpose() - Vstar).array()*inv_sigma)).rowwise().sum() ).array();
    U = fill_matrix(u,n);
    
    ustar = (K-U).colwise().maxCoeff();
    Ustar = fill_matrix(ustar,n);
    
    v = nu.array() + ustar.array() + sigma * Eigen::log( (Eigen::exp( (K - U - Ustar.transpose()).array()*inv_sigma)).transpose().rowwise().sum() ).array();
    V = fill_matrix(v,n);
    
    test = (u.array() - mu.array() - sigma * Eigen::log( (Eigen::exp( (K - V.transpose()).array()*inv_sigma)).rowwise().sum() ).array()) ;
    if(max_iter<iter_ipfp || test.rowwise().sum().abs().maxCoeff()<epsilon) {
      stop_condition = true;
    }
  }
  
  printf("iter IPFP: %lu \n", iter_ipfp);
  //std::cout << Ys << std::endl;
  
  Eigen::VectorXd A = Eigen::exp(-inv_sigma*u.array());
  Eigen::VectorXd B = Eigen::exp(-inv_sigma*v.array());
  
  double norm = A(0);
  double inv_norm = 1/norm;
  
  Eigen::VectorXd A_norm = inv_norm * A;
  Eigen::VectorXd B_norm = norm * B;
  
  Eigen::VectorXd a = -sigma * Eigen::log(A_norm.array());
  Eigen::VectorXd b = -sigma * Eigen::log(B_norm.array());
  
  return List::create(Named("a")=a, Named("b")=b, Named("iterIpfp")=iter_ipfp);
}
