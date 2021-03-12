#include "cgos.h"

/*
  Simply computes the values of cos(m*phi) for 0<=m<=n.
  Values are organised in a matrix.
  For 0<=i<phi.size() and 0<=j<=n:
  cp(i,j)=cos(j*phi(i))
*/

void cos_m_phi(MatDouble &cp,VecDouble &p)
{
  std::cout << "Computing the values of cos(m*phi)..." << std::endl;

  VecDouble m (cp.size2());
  std::copy(CountItDouble0,CountItDouble0+m.size(),m.begin());
  cp=boostublas::outer_prod(p,m);
  AutoMatOp(cp,[](double c) -> double { return cos(c); });
}
