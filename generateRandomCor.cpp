/*
 * Based on the matlab code for which I reproduce the heading:
 * 
 * References:
 %   [1] Mohsen Pourahmadi and Xiao Wang,
 %       Distribution of random correlation matrices: Hyperspherical parameterization of the Cholesky factor,
 %       Statistics & Probability Letters, Volume 106, November 2015, Pages 5-12
 %
 %   [2] Enes Makalic and Daniel F. Schmidt,
 %       An efficient algorithm for sampling from $\sin^k(x)$ for generating random correlation matrices
 %       arxiv, 2018
 % 
 %   
 * 
 * 
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace Eigen;

const float pi=3.1416;
namespace Eigen{

 namespace internal {
 template<typename Scalar>
 struct scalar_normal_dist_op
 {
   static std::mt19937 rng;                        // The uniform pseudo-random algorithm
   mutable std::normal_distribution<Scalar> norm; // gaussian combinator
  
   EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)
    
     template<typename Index>
     inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
     inline void seed(const uint64_t &s) { rng.seed(s); }
 };

 template<typename Scalar>
 std::mt19937 scalar_normal_dist_op<Scalar>::rng;

 template<typename Scalar>
 struct functor_traits<scalar_normal_dist_op<Scalar> >
 { enum { Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false }; };
 
 } // end namespace internal
}
VectorXf draw_theta(const VectorXf& k){
  unsigned int N; 
  float logconst;
  N = k.size();
  VectorXf x(N);
  x.setZero();
  Array<bool,Dynamic,1> accept;
  accept.setConstant(N,false);
  logconst = 2 * std::log(pi / 2);
  while( !accept.all()){
    for(int i=0; i < N ; i++){
      if(~accept(i)){
        float g1;
        float g2;
        g1 = R::rgamma(k(i)+1,1);
        g2 = R::rgamma(k(i)+1,1);
        x(i) = pi * (g1/(g1+g2));
      }
    }
    for(int i = 0 ; i < N ; i++){
      if(~accept(i))
        accept(i) = (std::log(R::runif(0,1)/k(i))  < logconst + std::log(std::sin(x(i))) - std::log(x(i)) - std::log(pi - x(i)) );
    }
  }
  return x;
}


//[[Rcpp::export]]
RowVectorXf cumprod(const RowVectorXf& X){
  int num_outer(X.size());
  RowVectorXf Y(X.size());
  float prod = 1;
  for(int i = 0;i<num_outer;i++)
  {
    prod *= std::sin(X(i));
    Y(i) = prod;
  }
  return Y;
}

// [[Rcpp::export]]
MatrixXf generateRandomCorL(unsigned int  p) {
  //Step 1 - generate angles theta from PDF (sin(theta))^k, k>=1, 0<theta<pi
  MatrixXf theta(p,p);
  double halfc;
  halfc = pi/180;
  theta.setZero();
  for (int j = 0 ; j < (p-1); j++)
    theta.col(j).segment( (j+1), (p-(j+1))) = draw_theta(  VectorXf(p-(j+1)).setConstant(p-(j+1)) );
  // Step 2 - construct lower triangular Cholesky factor
  MatrixXf L(p,p);
  L.setOnes();
  //L.triangularView<Eigen::Lower>().setConstant(1.0);
  MatrixXf tmp(p,p);
  tmp.setOnes();
  for (int i=1; i < p ; i++){
    L.row(i).segment(1,i)=cumprod(theta.row(i).segment(0, i));
    for(int j=0; j<i; j++ )
      tmp(i,j)=std::cos(theta(i,j)); //EIGEN COSINE FUNTION HAS ISSUES!!!!!
  }
  //return MatrixXf(p,p).setZero().selfadjointView<Eigen::Lower>().rankUpdate(MatrixXf(L.triangularView<Eigen::Lower>()).cwiseProduct(tmp) );
  return MatrixXf(L.triangularView<Eigen::Lower>()).cwiseProduct(tmp);
}

// [[Rcpp::export]]
MatrixXf generateRandomCor(unsigned int p){
  return MatrixXf(p,p).setZero().selfadjointView<Eigen::Lower>().rankUpdate(generateRandomCorL(p));
}

//[[Rcpp::export]]
MatrixXf generateMvn(unsigned int p,unsigned int N){
  MatrixXf L(p,p);
  Eigen::internal::scalar_normal_dist_op<float> randN;
  randN.seed(std::mt19937::default_seed);
  L = generateRandomCorL(p);
  MatrixXf X(N,p);
  X = L.triangularView<Eigen::Lower>() * Eigen::Matrix<float,Dynamic,-1>::NullaryExpr(p,N,randN);
  return X;
}

//[[Rcpp::export]]
MatrixXf generateMvnL(MatrixXf L, unsigned int N){
  Eigen::internal::scalar_normal_dist_op<float> randN;
  randN.seed(std::mt19937::default_seed);
  MatrixXf X(L.rows(),L.rows());
  X = L.triangularView<Eigen::Lower>() * Eigen::Matrix<float,Dynamic,-1>::NullaryExpr(L.rows(),N,randN);
  return X;
}

