#include "eit.h"

/*
  Computes  exp(-ix.\xi) at quadrature points x.
*/

void exp_minus_i_x_xi(VecComplex &expmixxi,const VecDouble &xi,MatDouble &x)
{
  VecDouble xix = -prec_inner_prod_vec_rows_mat(xi,x);
  expmixxi = complexfunvec(VecOp(xix,[](double c) -> double { return cos(c); }),VecOp(xix,[](double c) -> double { return sin(c); }));
}
