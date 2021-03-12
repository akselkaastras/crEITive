#include "pcc.h"

/*
  Computes the values of the densities at quadrature points,
  given their spherical harmonics series:
  f(q(z)) = sum f_n^m Y_n^m(z)
  where q is the map from the sphere to the surface.
*/

void dev_to_val(MatComplex &valdens,MatComplex &devdens,Conductivity &sigma)
{
  unsigned n=sigma.GetMaxDegree();
  MatDouble lf=sigma.GetAssociatedLegendreFunctions();
  MatComplex ep=sigma.GetComplexExponentials();
  unsigned ndponetwo=devdens.size2();

  unsigned ni,k;
  int l;

  valdens = ZeroMatComplex (2*(n+1)*(n+1),ndponetwo);

  for (unsigned ip=0;ip<2*(n+1);ip++) 
    {
      for (unsigned it=0;it<n+1;it++) 
	{
	  ni=it+(n+1)*ip;
	  for (unsigned kd=0;kd<ndponetwo;kd++)//degree and orders of solid harmonics (densities)
	    {
	      for (unsigned kl=0;kl<(n+1)*(n+1);kl++) //degree and orders of spherical harmonics (series)
		{
		  k=intsqrt(kl);//degree
		  l=kl-k*k-k;//order
		  valdens(ni,kd)+=devdens(kl,kd)*lf(it,(k*(k+1))/2+abs(l))*ep(ip,n+l);
		}
	    }
	}
    }
}
