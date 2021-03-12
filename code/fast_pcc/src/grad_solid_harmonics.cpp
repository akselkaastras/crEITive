#include "pcc.h"

/*
  Computes the gradient of the solid harmonics for non-negative
  degrees and non-negative orders at points x=(x1,x3,x3).
  They are organised in a matrix.
  For i>=0 and (k*(k+1))/2<=j<((k+1)*(k+2))/2, where 0<=k<=n
  R(i,j)=R_k^l(x(i)) where l=j-k*(k+1)/2

  The size of R is (x.size1(),((n+1)*(n+2))/2), where n is the maximal degree.

  We have the following formula for the solid harmonics:
  R_k^l(x) = r^k P_k^l(cos(theta)) exp(I*m*phi),
  Where (r,theta,phi) are the spherical coordinates of x and P_k^l are
  the associated Legendre functions
  P_k^l(t) = sqrt(((2*k+1)/(4*Pi))*fact(k-l)/fact(k+l))*(1-t^2)^(l/2) * d/dt^l P_k(t)
  Where P_k are the Legendre polynomials.

  /!\/!\ /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
  The leading Condon-Shortley phase term of (-1)^l is removed.
  /!\/!\ /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

  For l>=0 and k>=l+2
  R_k^l(x) = a_k^l * x3 * R_{k-1}^l(x)
  - b_k^l * (x1^2 + x2^2 + x3^2) * R_{k-2}^l(x)
  Where:
  a_k^l=sqrt((2*k+1)*(2*k-1)/((k+l)*(k-l)))
  b_k^l=sqrt((2*k+1)*(k+l-1)*(k-l-1)/((k+l)*(k-l)*(2*k-3)))

  dx1R_k^l(x) = a_k^l * x3 * dx1R_{k-1}^l(x)
  - b_k^l * ( 2 * x1 * R_{k-2}^l(x)
  + (x1^2 + x2^2 + x3^2) dx1R_{k-2}^l(x) )
  dx2R_k^l(x) = a_k^l * x3 * dx2R_{k-1}^l(x)
  - b_k^l * ( 2 * x2 * Y_{k-2}^l(x)
  + (x1^2 + x2^2 + x3^2) dx2R_{k-2}^l(x) )
  dx3R_k^l(x) = a_k^l * ( R_{k-1}^l(x) + x3 * dx3R_{k-1}^l(x) )
  - b_k^l * ( 2 * x3 * R_{k-2}^l(x)
  + (x1^2 + x2^2 + x3^2) dx3R_{k-2}^l(x) )

  The recusion formula is still valid for k>=1 and l=k-1
  if we define R_{k-1}^k by whatever we want for k>=0 since b_{k+1}^k=0.
  We then have R_k^{k-1}(x) = a_k^{k-1} * x3 * R_{k-1}^{k-1}(x) for k>=1.
  That is  R_k^{k-1}(x) = sqrt(2*k+1) * x3 * R_{k-1}^{k-1}(x).
  And
  dx1R_k^{k-1}(x) = sqrt(2*k+1) * x3 * dx1R_{k-1}^{k-1}(x)
  dx2R_k^{k-1}(x) = sqrt(2*k+1) * x3 * dx2R_{k-1}^{k-1}(x)
  dx3R_k^{k-1}(x) = sqrt(2*k+1) * ( R_{k-1}^{k-1}(x) + x3 * dx3R_{k-1}^{k-1}(x) )

  Seed values used for the recursion:
  R_0^0(x) = 1/sqrt(4*Pi)
  dx1R_0^0(x) = 0
  dx2R_0^0(x) = 0
  dx3R_0^0(x) = 0

  For l>=1:
  R_l^l(x) = sqrt((2*l+1)/(2*l)) * (x1 + I*x2) * R_{l-1}^{l-1}(x)
  dx1R_l^l(x) = sqrt((2*l+1)/(2*l)) * ( R_{l-1}^{l-1}(x) +
  (x1 + I*x2) * dx1R_{l-1}^{l-1}(x) )
  dx2R_l^l(x) = sqrt((2*l+1)/(2*l)) * ( I * R_{l-1}^{l-1}(x) +
  (x1 + I*x2) * dx2R_{l-1}^{l-1}(x) )
  dx3R_l^l(x) = sqrt((2*l+1)/(2*l)) * (x1 + I*x2) * dx3R_{l-1}^{l-1}(x)
*/

void grad_solid_harmonics(Array3Complex &gradR,MatDouble &x,unsigned &n)
{
  //std::cout << "Computing gradient of solid harmonics..." << std::endl;

  //For manipulation of beginning and ending columns for a degree
  unsigned colbegpp;
  unsigned colbegp;
  unsigned colendp;
  unsigned colbeg;
  unsigned colend;
  //Recursion coefficients
  double akl;
  double bkl;
  //Solid harmonics
  MatComplex R (x.size1(),((n+1)*(n+2))/2);
  //Misc
  double rx; //for x1^2+x2^2+x3^2
  dcomplex sx; //for x1+1i*x2

  for (unsigned i=0;i<x.size1();i++) //Fill the first three columns
    {
      R(i,0)=invroot4Pi;
      R(i,1)=sqrt(3.0)*invroot4Pi*x(i,2);
      R(i,2)=sqrt(3.0/2.0)*invroot4Pi*dcomplex(x(i,0),x(i,1));

      gradR[i][0][0]=0.0;
      gradR[i][0][1]=0.0;
      gradR[i][0][2]=0.0;

      gradR[i][1][0]=0.0;
      gradR[i][1][1]=0.0;
      gradR[i][1][2]=sqrt(3.0)*invroot4Pi;

      gradR[i][2][0]=sqrt(3.0/2.0)*invroot4Pi;
      gradR[i][2][1]=dcomplex(0.0,sqrt(3.0/2.0)*invroot4Pi);
      gradR[i][2][2]=0.0;
    }

  for (unsigned i=0;i<x.size1();i++) //Values of x
    {
      colbegpp=0;
      colbegp=1;
      colendp=2;
      colbeg=3;
      colend=5;

      rx=pow(x(i,0),2)+pow(x(i,1),2)+pow(x(i,2),2);
      sx=dcomplex(x(i,0),x(i,1));

      for (unsigned k=2;k<=n;k++) //Degrees of solid harmonics
	{
	  for (unsigned l=0;l<k-1;l++) //Order of solid harmonics: recursion formula
	    {
	      //Recursion coefficients
	      akl=sqrt((2.0*(double)k+1.0)*(2.0*(double)k-1.0)/(((double)k+(double)l)*((double)k-(double)l)));
	      bkl=sqrt((2.0*(double)k+1.0)*((double)k+(double)l-1.0)*((double)k-(double)l-1.0)/(((double)k+(double)l)*((double)k-(double)l)*(2.0*(double)k-3.0)));
	      //Recursion
	      R(i,colbeg+l)=akl*x(i,2)*R(i,colbegp+l)-bkl*rx*R(i,colbegpp+l);

	      gradR[i][colbeg+l][0]=akl*x(i,2)*gradR[i][colbegp+l][0]-bkl*(2.0*x(i,0)*R(i,colbegpp+l)+rx*gradR[i][colbegpp+l][0]);
	      gradR[i][colbeg+l][1]=akl*x(i,2)*gradR[i][colbegp+l][1]-bkl*(2.0*x(i,1)*R(i,colbegpp+l)+rx*gradR[i][colbegpp+l][1]);
	      gradR[i][colbeg+l][2]=akl*(R(i,colbegp+l)+x(i,2)*gradR[i][colbegp+l][2])-bkl*(2.0*x(i,2)*R(i,colbegpp+l)+rx*gradR[i][colbegpp+l][2]);

	    }
	  //Order k-1 and k
	  R(i,colend-1)=sqrt(2.0*(double)k+1.0)*x(i,2)*R(i,colendp);
	  R(i,colend)=sqrt((2.0*(double)k+1.0)/(2.0*(double)k))*sx*R(i,colendp);

	  gradR[i][colend-1][0]=sqrt(2.0*(double)k+1.0)*x(i,2)*gradR[i][colendp][0];
	  gradR[i][colend-1][1]=sqrt(2.0*(double)k+1.0)*x(i,2)*gradR[i][colendp][1];
	  gradR[i][colend-1][2]=sqrt(2.0*(double)k+1.0)*(R(i,colendp)+x(i,2)*gradR[i][colendp][2]);

	  gradR[i][colend][0]=sqrt((2.0*(double)k+1.0)/(2.0*(double)k))*(R(i,colendp)+sx*gradR[i][colendp][0]);
	  gradR[i][colend][1]=sqrt((2.0*(double)k+1.0)/(2.0*(double)k))*(Ic*R(i,colendp)+sx*gradR[i][colendp][1]);
	  gradR[i][colend][2]=sqrt((2.0*(double)k+1.0)/(2.0*(double)k))*sx*gradR[i][colendp][2];

	  //Change beginning and ending columns for the next degree
	  colbegpp=colbegp;
	  colbegp=colbeg;
	  colendp=colend;
	  colbeg=colend+1;
	  colend=colend+k+2;
	}
    }
}
