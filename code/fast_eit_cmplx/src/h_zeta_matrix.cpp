#include "eit.h"
#include "h_zeta.h"

/*
  Computes the matrix Mzeta of the Integral Operator with Kernel H_zeta on the unit sphere, where
  H_zeta(x) = G_zeta(x) - 1/(4*pi*|x|)
  G_zeta(x) being the Faddeev Green's function.

  In other words, Mzeta is the 2N^2x2N^2 matrix of f -> \int_{S^2} H_zeta(x-y) f(y) ds(y).

  The matrix Mzeta is computed for the quadrature points given by:
  x=(sqrt(1-t^2)*cos(phi),sqrt(1-t^2)*sin(phi),t).
  Where t are the n values given by the function gauss_legendre and phi=k*pi/N with 0<=k<=2*N-1.

  The quadrature points are organised by first increasing the index of t and then the one of phi.

  The matrix Mzeta is computed according to the quadrature rule (Wienert, Kress)
  \int_{S^2} H_zeta(x-y) f(y) ds(y) ~ \sum_{k=1}^N \sum_{l=0}^{2N-1} wtt(x_{kl}) f(x_{kl}) H_zeta(x-x_{kl})
  Where wtt is the weight at x_{kl}, it is given by wtt=\pi/N wt(theta_{_{kl}}) where wt are the weights given by the function gauss_legendre and where P_m are the Legendre polynomials.

  zeta = kappa*(k^T + i k)
  kappa >= 0
  k^T,k \in R^3 with k^T.k = 0 and |k|=|k^T|=1
*/

void h_zeta_matrix(MatDoubleColMaj &Mzeta,const double &kappa,const VecDouble &k,MatDouble &x,VecDouble &wtrep)
{
  unsigned n=wtrep.size();
  VecDouble ximxj (3);

  for (unsigned j=0;j<n;j++)
    {
      for (unsigned i=0;i<n;i++)
	{
	  ximxj=boostublas::row(x,i)-boostublas::row(x,j);
	  h_zeta(Mzeta(i,j),kappa,k,ximxj);
	  Mzeta(i,j)*=wtrep(j);
	}
    }
}
