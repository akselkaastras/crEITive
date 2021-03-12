#include "eit.h"

/*
  Computes the right hand side exp(ix.\zeta) at quadrature points x.

  zeta = kappa*(k^T + i k)
  kappa >= 0
  k^T,k \in R^3 with k^T.k = 0 and |k|=|k^T|=1

  Since the matrix of the boundary integral equation is real, the right hand side is stored as a two columns matrix, the first column for the real part and the second column for the imaginary part.
*/

void zeta_rhs(MatDoubleColMaj &rhs,const double &kappa,const VecDouble &k,const VecDouble &kT,MatDouble &x)
{
  VecDouble xk = VecOp(prec_inner_prod_vec_rows_mat(-kappa*k,x),[](double c) -> double { return exp(c); });
  VecDouble xkT = prec_inner_prod_vec_rows_mat(kappa*kT,x);

  boostublas::column(rhs,0)=boostublas::element_prod(VecOp(xkT,[](double c) -> double { return cos(c); }),xk);
  boostublas::column(rhs,1)=boostublas::element_prod(VecOp(xkT,[](double c) -> double { return sin(c); }),xk);
}
