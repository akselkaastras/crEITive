#include "radsym.h"

/*
  Simply computes the values of phi.
*/

void phi_values(VecDouble &p)
{
  std::cout << "Computing the values of phi..." << std::endl;

  unsigned n=p.size();

  for(unsigned k=0;k<n;k++) //Values of phi
    {
      p(k)=2.0*((double)k)*Pi/((double)n);
    }
}
