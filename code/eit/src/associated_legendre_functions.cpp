#include "eit.h"

/*
  Computes the values of the fully normalized associated Legendre functions
  for non-negative degrees and non-negative orders at points t.
  They are organised in a matrix.
  For i>=0 and (k*(k+1))/2<=j<((k+1)*(k+2))/2, where 0<=k<=n
  p(i,j)=P_k^l(t(i)) where l=j-k*(k+1)/2

  The size of p is (t.size(),((n+1)*(n+2))/2), where n is the maximal degree.

  We have the following formula for the fully normalized associated Legendre functions:
  P_k^l(t) = sqrt(((2*k+1)/(4*Pi))*fact(k-l)/fact(k+l))*(1-t^2)^(l/2) * d/dt^l P_k(t)
  Where P_k are the Legendre polynomials.

  /!\/!\ /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
  The leading Condon-Shortley phase term of (-1)^l is removed.
  /!\/!\ /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

  Boost library for the computations is not really appopriate since it computes non-normalized associated Legendre functions. I guess it gives better accuracy when using a recurrence formula on the normalized ones to compute the spherical harmonics instead of computing the non-normalized associated Legendre functions and dividing by the appropriate coefficient to get the spherical harmonics.
  Hence, this function is done by a three-term backward recursion formula for fixed l.

  For l>=0 and k>=l+2
  P_k^l(t) = a_k^l * t * P_{k-1}^l(t) - b_k^l * P_{k-2}^l(t)

  Where:
  a_k^l=sqrt((2*k+1)*(2*k-1)/((k+l)*(k-l)))
  b_k^l=sqrt((2*k+1)*(k+l-1)*(k-l-1)/((k+l)*(k-l)*(2*k-3)))

  The recusion formula is still valid for k>=1 and l=k-1
  if we define P_{k-1}^k by whatever we want for k>=0 since b_{k+1}^k=0.
  We then have P_k^{k-1}(t) = a_k^{k-1} * t * P_{k-1}^{k-1}(t) for k>=1.
  That is  P_k^{k-1}(t) = sqrt(2*k+1) * t * P_{k-1}^{k-1}(t).

  Seed values used for the recursion:
  P_0^0(t) = 1/sqrt(4*Pi)
  P_l^l(t) = sqrt((2*l+1)/(2*l)) * sqrt(1-t^2) * P_{l-1}^{l-1}(t) for l>=1.
*/

void associated_legendre_functions(MatDouble &lf,VecDouble &t)
{
  std::cout << "Computing fully normalized associated Legendre functions..." << std::endl;

  unsigned n=t.size();
  //For manipulation of beginning and ending columns for a degree
  unsigned colbegpp;
  unsigned colbegp;
  unsigned colendp;
  unsigned colbeg;
  unsigned colend;
  //Recursion coefficients
  double akl;
  double bkl;

  for (unsigned i=0;i<n;i++) //Fill the first three columns
    {
      lf(i,0)=invroot4Pi;
      lf(i,1)=sqrt(3.0)*t(i)*lf(i,0);
      lf(i,2)=sqrt(3.0/2.0)*sqrt(1.0-pow(t(i),2.0))*lf(i,0);
    }

  for (unsigned i=0;i<n;i++) //Values of t
    {
      colbegpp=0;
      colbegp=1;
      colendp=2;
      colbeg=3;
      colend=5;

      for (unsigned k=2;k<n;k++) //Degrees of Legendre functions
	{
	  for (unsigned l=0;l<k-1;l++) //Order of Legendre functions: recursion formula
	    {
	      //Recursion coefficients
	      akl=sqrt((2.0*(double)k+1.0)*(2.0*(double)k-1.0)/(((double)k+(double)l)*((double)k-(double)l)));
	      bkl=sqrt((2.0*(double)k+1.0)*((double)k+(double)l-1.0)*((double)k-(double)l-1.0)/(((double)k+(double)l)*((double)k-(double)l)*(2.0*(double)k-3.0)));
	      //Recursion
	      lf(i,colbeg+l)=akl*t(i)*lf(i,colbegp+l)-bkl*lf(i,colbegpp+l);
	    }
	  //Order k-1 and k
	  lf(i,colend-1)=sqrt(2.0*(double)k+1.0)*t(i)*lf(i,colendp);
	  lf(i,colend)=sqrt((2.0*(double)k+1.0)/(2.0*(double)k))*sqrt(1.0-pow(t(i),2.0))*lf(i,colendp);
	  //Change beginning and ending columns for the next degree
	  colbegpp=colbegp;
	  colbegp=colbeg;
	  colendp=colend;
	  colbeg=colend+1;
	  colend=colend+k+2;
	}
    }
}
