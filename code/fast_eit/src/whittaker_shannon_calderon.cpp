#include "eit.h"

/*
  Computes the conductivity by the Calderón approximation at points x using the Whittaker-Shannon interpolation formula:
  sigma(x) = 1 - 1/(4*Pi^2) \sum qhat(ii*Pi,jj*Pi,kk*Pi)/(ii^2+jj^2+kk^2) exp(Ic*(ii*x(0)+jj*x(1)+kk*x(2))*Pi) 
*/

void whittaker_shannon_calderon(VecDouble &sig,MatDouble &x,MatDouble &xi,Array3Complex &qhat)
{
  std::cout << "Computing the conductivity..." << std::endl;

  unsigned ng3 = xi.size1();
  unsigned ng = intsqrt3(ng3);
  VecDouble normxi2 = vec_norm_2_square_rows_mat(xi);

  //Avoid dividing by zero
  normxi2((ng3-1)/2)=1.0;
  //Put qhat(0,0,0)=0
  qhat[(ng-1)/2][(ng-1)/2][(ng-1)/2]=0.0;

  //Compute sigma
  VecComplex qhatline (ng3);
  std::copy(qhat.origin(),qhat.origin()+qhat.num_elements(),qhatline.begin());

  MatDouble xix = prec_inner_prod_mat_mat_rows(xi,x);
  sig = -real(boostublas::prec_prod(boostublas::element_div(qhatline,normxi2),complexfunmat(MatOp(xix,[](double c) -> double { return cos(c); }),MatOp(xix,[](double c) -> double { return sin(c); }))))/(4.0*square(Pi));
  AutoVecOp(sig,plus_one<double>);
}
