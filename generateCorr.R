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

