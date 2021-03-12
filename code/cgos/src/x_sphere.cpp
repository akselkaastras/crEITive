#include "cgos.h"

/*
  Computes the quadrature points on the sphere in cartesian coordinates.
  x=(sqrt(1-t^2)*cos(phi),sqrt(1-t^2)*sin(phi),t)
  Where t are the N values given by the function gauss_legendre and phi=k*pi/N with 0<=k<=2*N-1.

  The quadrature points are organised by first increasing the index of t and then the one of phi.
*/

void x_sphere(MatDouble &x,VecDouble &t,VecDouble &phi)
{
  unsigned n=t.size();

  VecDouble sit = VecOp((VecDouble)(ScalVecDouble(n,1.0)-VecOp(t,square<double>)),[](double c) -> double { return sqrt(c); });

  boostublas::column(x,0) = kron_prod_vec(VecOp(phi,[](double c) -> double { return cos(c); }),sit);
  boostublas::column(x,1) = kron_prod_vec(VecOp(phi,[](double c) -> double { return sin(c); }),sit);
  boostublas::column(x,2) = rep_vec_1(t,2*n);
}
