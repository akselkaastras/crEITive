#include "radsym.h"

/*
  Computes q on a radius.
*/

void q_on_radius(VecDouble &r,VecDouble &valq,MatDouble &sigma)
{
  unsigned nsigma=sigma.size1();
  unsigned n=r.size();
  valq = ZeroVecDouble (n);
  VecDouble valsigma = ZeroVecDouble (n); //sigma-1 which will become sigma
  VecDouble valdsigma = ZeroVecDouble (n);
  VecDouble vald2sigma = ZeroVecDouble (n);
  double h=1.0/((double)n-1.0);
  double ri,re,amp,var,var2,fth,fth2;
  double valfunc,valfuncp,valfuncpp;
  unsigned mmin,mmax;

  //Discretization points
  for (unsigned m=0;m<n;m++)
    {
      r(m)=(double)m*h;
    }

  //Compute values of q
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
	  fth2=pow(fth,2);
	  valfunc=(amp-1.0)/2.0*(1.0-fth);
	  valfuncp=2.0*(1.0-amp)/(re-ri)*(1.0+var2)/pow(1.0-var2,2)*(1.0-fth2);
	  valfuncpp=8.0*(1.0-amp)/pow(re-ri,2)*(var*(3.0+var2)*(1.0-var2)-2.0*pow(1.0+var2,2)*fth)/pow(1.0-var2,4)*(1.0-fth2);

	  valsigma(m)+=valfunc;
	  valdsigma(m)+=valfuncp;
	  vald2sigma(m)+=valfuncpp+2.0/r(m)*valfuncp;
	}
    }
  for (unsigned m=0;m<n;m++)
    {
      valsigma(m)+=1.0;
      valq(m)=vald2sigma(m)/(2.0*valsigma(m))-pow(valdsigma(m),2)/(4.0*pow(valsigma(m),2));
    }
}
