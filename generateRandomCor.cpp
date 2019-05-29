#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace Eigen;

const float pi=3.1416;
VectorXf draw_theta(const VectorXf& k){
  unsigned int N; 
  float logconst;
  N = k.size();
  VectorXf x(N);
  x.setZero();
  Array<bool,Dynamic,1> accept;
  accept.setConstant(N,false);
  logconst = 2 * log(pi / 2);
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
        accept(i) = (log(R::runif(0,1)/k(i))  < logconst + log(sin(x(i))) - log(x(i)) - log(pi - x(i)) );
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
      prod *= sin(X(i));
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
      tmp(i,j)=cos(theta(i,j)); //EIGEN COSINE FUNTION HAS ISSUES!!!!!
  }
  //return MatrixXf(p,p).setZero().selfadjointView<Eigen::Lower>().rankUpdate(MatrixXf(L.triangularView<Eigen::Lower>()).cwiseProduct(tmp) );
  return MatrixXf(L.triangularView<Eigen::Lower>()).cwiseProduct(tmp);
}

// [[Rcpp::export]]
MatrixXf generateRandomCor(unsigned int p){
  return MatrixXf(p,p).setZero().selfadjointView<Eigen::Lower>().rankUpdate(generateRandomCorL(p));
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
set.seed(123)
#D<-generateRandomCor(5)
Lc <- generateRandomCor(5)

# lets try to generate  then theta matrix ourselves

generate_theta <- function( k){
  N= length(k)
  logconst = 2 * log(pi/ 2);
  accept = rep(FALSE,N)
  x = rep(NA,N);
  while( !all(accept)){
    iX = !accept
    g1 =  rgamma(sum(iX),shape = k[iX]+1,scale=1)
    g2 = rgamma(sum(iX),shape = k[iX]+1,scale=1)
      
    x[iX] = pi *(g1/(g1+g2))
   
    accept[iX] = log(runif(n = sum(iX)))/k[iX]  < logconst + log(sin(x[iX])) - log(x[iX]) - log(pi - x[iX]) ;
  }
  x
}
theta <- matrix(rep(0,5^2),ncol=5) 
#theta <-D 
p=5
e <- rep(1,5)
for (j in 1:(p-1)){
#  theta[(j+1):p,j] = generate_theta( (p-j)*e[(j+1):p] );
}
L = matrix(rep(1,p^2),ncol=p)
for( i in 2:p){
   L[i,2:i] = base::cumprod( sin(theta[i,1:i-1]) );
}
cosL <- cos(theta)
cosL[upper.tri(cosL)] <-0
L <- L*cosL

C = L%*%t(L);
print(C)

  */
