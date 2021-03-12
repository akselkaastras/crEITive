#include "pcc.h"
#include "associated_legendre_functions.h"
#include "exp_i_m_phi.h"

/*
  Computes the matrix of the operator K':
  (K'\varphi)(z) = \int_{surface} \frac{\partial\Phi(z,z')}{\partial \nu(z)}
  See 3.6 of "Inverse Acoustic and Electromagnetic Scattering Theory" by
  David Colton and Rainer Kress.
*/

void kp_matrix(MatComplex &Kp,Conductivity &sigma)
{
  unsigned n=sigma.GetMaxDegree();
  MatDouble lf=sigma.GetAssociatedLegendreFunctions();
  MatComplex ep=sigma.GetComplexExponentials();
  MatDouble cttxmoy=sigma.GetCosElevationRotationsSpherePoints();
  MatDouble ptxmoy=sigma.GetAzimuthRotationsSpherePoints();
  MatDouble kernf=sigma.GetKernelF();
  std::cout << "Test: Is Kernel nnan? ...... " << kernf(1,1) << std::endl;
  unsigned nponetwo=pow(n+1,2);
  unsigned ni,nj,k,kprime;
  int l,lprime;

  MatDouble lftxmoy (2*(n+1)*(n+1),((n+1)*(n+2))/2);
  MatComplex eptxmoy (2*(n+1)*(n+1),2*n+1);

  Kp = ZeroMatComplex (Kp.size1(),Kp.size2());

  VecDouble cttxmoyrow (2*(n+1)*(n+1));
  VecDouble ptxmoyrow (2*(n+1)*(n+1));

  //compute Kp
  for (ni=0;ni<2*(n+1)*(n+1);ni++) 
    {
      cttxmoyrow=boostublas::row(cttxmoy,ni);
      ptxmoyrow=boostublas::row(ptxmoy,ni);
      associated_legendre_functions(lftxmoy,cttxmoyrow,n);
      exp_i_m_phi(eptxmoy,ptxmoyrow);

#pragma omp parallel default(none) shared(Kp,kernf,lf,lftxmoy,ep,eptxmoy,n,nponetwo,ni) private(k,l,kprime,lprime,nj)
      {
#pragma omp for collapse(2)
	for (unsigned kl=0;kl<nponetwo;kl++) 
	  {
	    for (unsigned klprime=0;klprime<nponetwo;klprime++) 
	      {
		k=intsqrt(kl);//degree
		l=kl-k*k-k;//order
		kprime=intsqrt(klprime);//degree
		lprime=klprime-kprime*kprime-kprime;//order

		//make the sum
		for (unsigned jp=0;jp<2*(n+1);jp++) 
		  {
		    for (unsigned jt=0;jt<n+1;jt++) 
		      {
			nj=jt+(n+1)*jp;
			Kp(kl,klprime)+=kernf(ni,nj)*lf(jt,(k*(k+1))/2+abs(l))*lftxmoy(nj,(kprime*(kprime+1))/2+abs(lprime))*ep(jp,n-l)*eptxmoy(nj,n+lprime);
      

		      }
		  }

	      }
	  }
      }
    }
  Kp/=4.0*Pi;
}
