#include "radsym.h"

/*
  Simply computes the values of cos(m*phi) for 0<=m<=n.
  Values are organised in a matrix.
  For 0<=i<phi.size() and 0<=j<=n:
  cp(i,j)=cos(j*phi(i))
*/

void cos_m_phi(MatDouble &cp,VecDouble &p)
{
  std::cout << "Computing the values of cos(m*phi)..." << std::endl;

  for (unsigned m=0;m<cp.size2();m++)
    {
      for (unsigned k=0;k<cp.size1();k++) //Values of phi
	{
	  cp(k,m)=cos((double)m*p(k));
	}
    }
}
