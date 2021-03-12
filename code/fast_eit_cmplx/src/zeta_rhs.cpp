#include "eit.h"

/*
  Computes the right hand side exp(ix.\zeta) at quadrature points x.

  zeta = kappa*(k^T + i k)
  kappa >= 0
  k^T,k \in R^3 with k^T.k = 0 and |k|=|k^T|=1
*/

void zeta_rhs(VecComplex &rhs,const double &kappa,const VecDouble &k,const VecDouble &kT,MatDouble &x)
{
  VecDouble xk = VecOp(prec_inner_prod_vec_rows_mat(-kappa*k,x),[](double c) -> double { return exp(c); });
  VecDouble xkT = prec_inner_prod_vec_rows_mat(kappa*kT,x);
  rhs=complexfunvec(boostublas::element_prod(VecOp(xkT,[](double c) -> double { return cos(c); }),xk),boostublas::element_prod(VecOp(xkT,[](double c) -> double { return sin(c); }),xk));
}
