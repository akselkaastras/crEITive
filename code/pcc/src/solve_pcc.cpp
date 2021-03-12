#include "pcc.h"

/*
  Solves the boundary intgeral equation. No use of Lapack since the matrix does not seem to have any special property.
*/

void solve_pcc(MatComplex &KpKbarp,MatComplex &rhs)
{
  std::cout << "Solving linear system..." << std::endl;

  unsigned n=KpKbarp.size1();

  // LU factorization
  PermMat pm (n);
  int info=lu_factorize(KpKbarp,pm);

  //Solve the boundary integral equation
  swap_rows(pm,rhs); //Permute rows of right hand side
  rhs=boostublas::solve(KpKbarp,boostublas::solve(KpKbarp,rhs,boostublas::unit_lower_tag()),boostublas::upper_tag()); //Solve system
}
