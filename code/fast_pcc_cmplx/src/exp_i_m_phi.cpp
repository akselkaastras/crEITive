#include "pcc.h"

/*
  Simply computes the values of exp(Ic*m*phi) for -n<=m<=n.
  Values are organised in a matrix.
  For 0<=i<phi.size() and 0<=j<=2*n:
  ep(i,j)=exp(Ic*(j-n)*phi(i))
*/

void exp_i_m_phi(MatComplex &ep,VecDouble &p)
{
  //std::cout << "Computing the values of exp(Ic*m*phi)..." << std::endl;

  unsigned n = (ep.size2()-1)/2;

  for (unsigned k=0;k<ep.size1();k++) //Values of phi
    {
      ep(k,n)=dcomplex(1.0,0.0);
      for (unsigned m=1;m<=n;m++)
	{
	  ep(k,n+m)=exp(dcomplex(0.0,(double)m*p(k)));
	  ep(k,n-m)=conj(ep(k,n+m));
	}
    }
}
