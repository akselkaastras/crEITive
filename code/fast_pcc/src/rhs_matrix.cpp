#include "pcc.h"
#include "grad_solid_harmonics.h"

/*
  Computes the matrix of the right hand side (without the factor 2 mu):
  \int_S (\partial_\nu R_{n'}^{m'}) (q(z)) Y_n^{-m} (z) ds(z)
  Where R_{n'}^{m'} are the solid harmonics, Y_n^{-m} the spherical harmonics,
  q is the map from the sphere S to the surface.
*/

void rhs_matrix(MatComplex &rhs,Conductivity &sigma,unsigned &nd)
{

  unsigned n=sigma.GetMaxDegree();
  double width=sigma.GetWidth();
  if (width < 0.01)
  {
  VecDouble alpha=sigma.GetAlpha();
  MatDouble x=sigma.GetSurfacePoints();
  MatDouble nu=sigma.GetOutwardNormal();
  MatDouble lf=sigma.GetAssociatedLegendreFunctions();
  MatComplex ep=sigma.GetComplexExponentials();

  Array3Complex gradRd (boost::extents[2*(n+1)*(n+1)][((nd+1)*(nd+2))/2][3]);

  unsigned nponetwo=(n+1)*(n+1);
  unsigned ndponetwo=(nd+1)*(nd+1);

  unsigned ni,k,kd;
  int l,ld;

  grad_solid_harmonics(gradRd,x,nd);

  rhs = ZeroMatComplex (rhs.size1(),rhs.size2());

#pragma omp parallel default(none) shared(rhs,gradRd,alpha,lf,ep,nu,nponetwo,ndponetwo,n) private(k,l,kd,ld,ni)
  {
#pragma omp for collapse(2)
    for (unsigned kl=0;kl<nponetwo;kl++) 
      {
	for (unsigned kld=0;kld<ndponetwo;kld++) 
	  {
	    k=intsqrt(kl);//degree
	    l=kl-k*k-k;//order
	    kd=intsqrt(kld);//degree
	    ld=kld-kd*kd-kd;//order
	    VecComplex gradRdi (3);
	    dcomplex dnuRdi;

	    //make the sum
	    for (unsigned ip=0;ip<2*(n+1);ip++) 
	      {
		for (unsigned it=0;it<n+1;it++) 
		  {
		    ni=it+(n+1)*ip;
		    gradRdi(0)=gradRd[ni][(kd*(kd+1))/2+abs(ld)][0];
		    gradRdi(1)=gradRd[ni][(kd*(kd+1))/2+abs(ld)][1];
		    gradRdi(2)=gradRd[ni][(kd*(kd+1))/2+abs(ld)][2];
		    if (ld<0)
		      {
			gradRdi=conj(gradRdi);
		      }
		    dnuRdi=boostublas::prec_inner_prod(gradRdi,boostublas::row(nu,ni));
		    rhs(kl,kld)+=alpha(it)*lf(it,(k*(k+1))/2+abs(l))*ep(ip,n-l)*dnuRdi;
		  }

	      }
	  }
      }
  }
}
 else
 {

	  VecDouble alpha=sigma.GetAlpha();
	  MatDouble x=sigma.GetSurfacePoints();
	  MatDouble x2=sigma.GetSurfacePoints2();
	  MatDouble nu=sigma.GetOutwardNormal();
	  MatDouble nu2=sigma.GetOutwardNormal2();
	  MatDouble lf=sigma.GetAssociatedLegendreFunctions();
	  MatComplex ep=sigma.GetComplexExponentials();

	  Array3Complex gradRd (boost::extents[2*(n+1)*(n+1)][((nd+1)*(nd+2))/2][3]);
	  Array3Complex gradRd2 (boost::extents[2*(n+1)*(n+1)][((nd+1)*(nd+2))/2][3]);

	  unsigned nponetwo=(n+1)*(n+1);
	  unsigned ndponetwo=(nd+1)*(nd+1);

	  unsigned ni,k,kd;
	  int l,ld;

	  grad_solid_harmonics(gradRd,x,nd);
	  grad_solid_harmonics(gradRd2,x2,nd);

	  rhs = ZeroMatComplex (rhs.size1(),rhs.size2());

	#pragma omp parallel default(none) shared(rhs,gradRd,gradRd2,alpha,lf,ep,nu,nu2,nponetwo,ndponetwo,n) private(k,l,kd,ld,ni)
	  {
	#pragma omp for collapse(2)
	    for (unsigned kl=0;kl<nponetwo;kl++) 
	      {
		for (unsigned kld=0;kld<ndponetwo;kld++) 
		  {
		    k=intsqrt(kl);//degree
		    l=kl-k*k-k;//order
		    kd=intsqrt(kld);//degree
		    ld=kld-kd*kd-kd;//order
		    VecComplex gradRdi (3);
		    VecComplex gradRdi2 (3);
		    dcomplex dnuRdi;
		    dcomplex dnuRdi2;

		    //make the sum
		    for (unsigned ip=0;ip<2*(n+1);ip++) 
		      {
			for (unsigned it=0;it<n+1;it++) 
			  {
			    ni=it+(n+1)*ip;
			    gradRdi(0)=gradRd[ni][(kd*(kd+1))/2+abs(ld)][0];
			    gradRdi(1)=gradRd[ni][(kd*(kd+1))/2+abs(ld)][1];
			    gradRdi(2)=gradRd[ni][(kd*(kd+1))/2+abs(ld)][2];
			    gradRdi2(0)=gradRd2[ni][(kd*(kd+1))/2+abs(ld)][0];
			    gradRdi2(1)=gradRd2[ni][(kd*(kd+1))/2+abs(ld)][1];
			    gradRdi2(2)=gradRd2[ni][(kd*(kd+1))/2+abs(ld)][2];
			    if (ld<0)
			      {
				gradRdi=conj(gradRdi);
				gradRdi2=conj(gradRdi2);
			      }
			    dnuRdi=boostublas::prec_inner_prod(gradRdi,boostublas::row(nu,ni));
			    dnuRdi2=boostublas::prec_inner_prod(gradRdi2,boostublas::row(nu2,ni));
			    rhs(kl,kld)+=alpha(it)*lf(it,(k*(k+1))/2+abs(l))*ep(ip,n-l)*dnuRdi;
			    rhs(kl,kld)+=alpha(it)*lf(it,(k*(k+1))/2+abs(l))*ep(ip,n-l)*dnuRdi2;
			  }

		      }
		  }
	      }
	  }

 }



}
