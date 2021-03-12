#include "radsym.h"

/*
  Computes the conductivity on a radius.
*/

void conductivity_on_radius(VecDouble &r,VecDouble &valsigma,MatDouble &sigma)
{
  unsigned nsigma=sigma.size1();
  unsigned n=r.size();
  valsigma = ZeroVecDouble (n);
  double h=1.0/((double)n-1.0);
  double ri,re,amp,var,var2,fth;
  unsigned mmin,mmax;

  //Discretization points
  for (unsigned m=0;m<n;m++)
    {
      r(m)=(double)m*h;
    }

  //Compute values of sigma-1
  for (unsigned k=0;k<nsigma;k++)
    {
      re=sigma(k,0);
      ri=sigma(k,0)*sigma(k,2);
      amp=sigma(k,1);
      mmin=floor(ri/h)+1;
      mmax=ceil(re/h);
      for (unsigned m=0;m<mmin;m++)
	{
	  valsigma(m)+=amp-1.0;
	}
      for (unsigned m=mmin;m<mmax;m++)
	{
	  var=2.0*(r(m)-ri)/(re-ri)-1.0;
	  var2=pow(var,2);
	  fth=tanh(2.0*var/(1.0-var2));
	  valsigma(m)+=(amp-1.0)/2.0*(1.0-fth);
	}
      //zero for m>=mmax, done by initialization of valsigma
    }

  // Add 1
  for (unsigned m=0;m<n;m++)
    {
      valsigma(m)+=1.0;
    }
}
