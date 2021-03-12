#include "eit.h"

/*
  Computes q at points x using the Whittaker-Shannon interpolation formula:
  q(x) = 1/8 \sum qhat(ii*Pi,jj*Pi,kk*Pi) exp(Ic*(ii*x(0)+jj*x(1)+kk*x(2))*Pi) 
*/

void whittaker_shannon(VecDouble &q,MatDouble &x,MatDouble &xi,Array3Complex &qhat)
{
  std::cout << "Computing q..." << std::endl;

  unsigned ng3=xi.size1();

  VecComplex qhatline (ng3);
  std::copy(qhat.origin(),qhat.origin()+qhat.num_elements(),qhatline.begin());

  MatDouble xix = prec_inner_prod_mat_mat_rows(xi,x);
  q = real(boostublas::prec_prod(qhatline,complexfunmat(MatOp(xix,[](double c) -> double { return cos(c); }),MatOp(xix,[](double c) -> double { return sin(c); }))))/8.0;
}
