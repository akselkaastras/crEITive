#include "eit.h"

/*
  Simply computes the values of phi.
*/

void phi_values(VecDouble &p)
{
  std::cout << "Computing the values of phi..." << std::endl;

  std::copy(CountItDouble0,CountItDouble0+p.size(),p.begin());
  p*=2.0*Pi/((double)p.size());
}
