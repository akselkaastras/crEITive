#include "eit.h"

/*
  Computes the scattering transform
  t(\xi,\zeta)=\int_{S^2} exp(-i*x.(\xi+\zeta) [(\Lambda_\sigma-\Lambda_0)\psi_\zeta](x) ds(x)

  The quadrature points are given by x=(sqrt(1-t^2)*cos(phi),sqrt(1-t^2)*sin(phi),t).
  Where t are the n values given by the function gauss_legendre and phi=k*pi/N with 0<=k<=2*N-1.

  The quadrature points are organised by first increasing the index of t and then the one of phi.

  zeta = kappa*(k^T + i k)
  kappa >= 0
  k^T,k \in R^3 with k^T.k = 0 and |k|=|k^T|=1
*/

void scattering_transform(dcomplex &txizeta,VecComplex &expmixxipzeta,MatDoubleColMaj &diffdtn,MatDoubleColMaj &psizeta,VecDouble &wtrep)
{
  // txizeta=boostublas::sum(boostublas::element_prod(wtrep,boostublas::element_prod(expmixxipzeta,complexfunmat2((MatDouble)boostublas::prec_prod(diffdtn,psizeta)))));
  // Blas product (faster)
  char trans[] = "N";
  double alpha=1.0;
  double beta=0.0;
  int n=wtrep.size();
  int one=1;
  MatDoubleColMaj pdp (n,2);
  dgemv_(trans,n,n,alpha,&diffdtn(0,0),n,&psizeta(0,0),one,beta,&pdp(0,0),one);
  dgemv_(trans,n,n,alpha,&diffdtn(0,0),n,&psizeta(0,1),one,beta,&pdp(0,1),one);
  txizeta=boostublas::sum(boostublas::element_prod(wtrep,boostublas::element_prod(expmixxipzeta,complexfunmat2(pdp))));
}