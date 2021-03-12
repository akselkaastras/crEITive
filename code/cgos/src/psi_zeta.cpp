#include "cgos.h"

/*
  Computes the solution \psi_\zeta of the boundary integral equation:
  [I+S_\zeta(\Lambda_\sigma-\Lambda_0)]\psi_\zeta=exp(i*x\zeta)

  Fzeta will denote I+S_\zeta(\Lambda_\sigma-\Lambda_0).

  The quadrature points are given by x=(sqrt(1-t^2)*cos(phi),sqrt(1-t^2)*sin(phi),t).
  Where t are the n values given by the function gauss_legendre and phi=k*pi/N with 0<=k<=2*N-1.

  The quadrature points are organised by first increasing the index of t and then the one of phi.

  zeta = kappa*(k^T + i k)
  kappa >= 0
  k^T,k \in R^3 with k^T.k = 0 and |k|=|k^T|=1

  On input, psizeta is the right hand side, on output, it is the solution.
*/

void psi_zeta(MatDoubleColMaj &Fzeta,MatDoubleColMaj &psizeta)
{
  //Lapack
  int n = psizeta.size1(); //Size of the system
  int nrhs = 2; //Number of right hand sides
  int lda = n; //Leading dimension of Fzeta
  int ldb = n; //Leading dimension of psizeta
  VecInt ipiv (n); //Pivot indices
  int info; //Information integer
  dgesv_(n,nrhs,&Fzeta(0,0),lda,&ipiv(0),&psizeta(0,0),ldb,info);
}
