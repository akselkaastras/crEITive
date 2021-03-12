#include "pcc.h"

/*
  Computes the normal derivative of the single layer potential.
  - \int_surface (1-|y|^2)/(4*Pi*|x-y|^3) f(y) ds(y)
  where f is the density.
  This is computed for the quadrature points x on the sphere and the densities f
  computed for each solid harmonic.
*/

void normal_derivatives(MatComplex &ND,MatComplex &valdens,VecConductivity &sigma,MatDouble &xd)
{
  std::cout << "Computing the normal derivatives..." << std::endl;

  unsigned nsig=sigma.size();
  unsigned twondponetwo=xd.size1();
  unsigned ndponetwo=twondponetwo/2;

  unsigned kstart,kend,ni;

  ND = ZeroMatComplex (twondponetwo,ndponetwo);

#pragma omp parallel default(none) shared(ND,sigma,xd,valdens,twondponetwo,ndponetwo,nsig) private(kstart,kend,ni)
  {
#pragma omp for collapse(2)
    for (unsigned kx=0;kx<twondponetwo;kx++)//quadrature points on the sphere
      {
	for (unsigned kd=0;kd<ndponetwo;kd++)//degree and orders of spherical harmonics
	  {
	    VecDouble xk = boostublas::row(xd,kx);
	    VecComplex valdenskd=boostublas::column(valdens,kd);

	    kstart=0;
	    for (unsigned k=0;k<nsig;k++)//conductivities
	      {
		unsigned n=sigma(k).GetMaxDegree();
		kend=kstart+2*(n+1)*(n+1);
		VecDouble alpha=sigma(k).GetAlpha();
		MatDouble x=sigma(k).GetSurfacePoints();
		VecDouble jx=sigma(k).GetJacobian();
		VecComplex valdensk=VecRangeComplex(valdenskd,boostublas::range(kstart,kend));

		//make the sum (quadrature points)
		for (unsigned ip=0;ip<2*(n+1);ip++) 
		  {
		    for (unsigned it=0;it<n+1;it++) 
		      {
			ni=it+(n+1)*ip;
			VecDouble xi=boostublas::row(x,ni);
			double normxi2=pow(xi(0),2)+pow(xi(1),2)+pow(xi(2),2);
			double normdiffx=boostublas::norm_2(xk-xi);
			ND(kx,kd)+=alpha(it)*jx(ni)*valdensk(ni)*(normxi2-1.0)/pow(normdiffx,3);
		      }
		  }
		kstart=kend;
	      }
	  }
      }
  }
  ND/=4.0*Pi;
}
