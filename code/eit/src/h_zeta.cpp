#include "eit.h"

/*
  Computes the harmonic function H_zeta(x) = G_zeta(x) - 1/(4*pi*|x|)
  Where G_zeta(x) is the Faddeev Green's function. It is computed using the following formula.

  G_(zeta)(x) = exp(-kappa*r)/(4*pi*r) - kappa/(4*pi) \int_{k.hatx}^1 exp(-kappa*r*u)/(sqrt(1-u^2))*J_1(kappa*r*sqrt(1-u^2)) du

  zeta = kappa*(k^T + i k)
  kappa >= 0
  k^T,k \in R^3 with k^T.k = 0 and |k|=|k^T|=1
  r = |x|
  hatx = x/|x|
  J_1 is the Bessel function of the first kind and of order 1.

  H_(zeta)(x) = - kappa/(4*pi)*[exprel(-kappa*r)+\int_{k.hatx}^1 exp(-kappa*r*u)/(sqrt(1-u^2))*J_1(kappa*r*sqrt(1-u^2)) du]

  exprel(u)=(exp(u)-1)/u

  Note: J_1(kappa*r*sqrt(1-u^2))/(sqrt(1-u^2)) is replaced by
  kappa*r*(J_0(kappa*r*sqrt(1-u^2))+J_2(kappa*r*sqrt(1-u^2)))/2
  to avoid having to consider u=-1 or u=1.
*/

double kern_h_zeta(double u, void * params)
{
  double kappar = *(double *) params;
  double root1mu2=sqrt(1.0-u*u);
  double kern_h_zeta;
  // if (u+1.0<=eps) //u equal to -1
  //   {
  //     kern_h_zeta=0.5*kappar*exp(kappar);
  //   }
  // else if (1.0-u<=eps) //u equal to 1
  //   {
  //     kern_h_zeta=0.5*kappar*exp(-kappar);
  //   }
  // else
  //   {
  //     kern_h_zeta=exp(-kappar*u)/root1mu2*boostmath::cyl_bessel_j(1,kappar*root1mu2);
  //   }
  kern_h_zeta=0.5*kappar*exp(-kappar*u)*(boostmath::cyl_bessel_j(0,kappar*root1mu2)+boostmath::cyl_bessel_j(2,kappar*root1mu2));
  return kern_h_zeta;
}

void h_zeta(double &hzeta,const double &kappa,const VecDouble &k,VecDouble &x)
{
  double r=boostublas::norm_2(x); //||x||
  double kappar=kappa*r;
  VecDouble hatx; //x/||x||
  double s; //k.hatx

  double epsabs=0.0;
  double epsrel=sqrt(eps);
  size_t neval;
  double error;
  gsl_function F;

  // Disable error, if relative error is not matched, last value is taken,
  // that is, the one computed with 87 points
  gsl_set_error_handler_off ();

  F.function = &kern_h_zeta;
  F.params = &kappar;

  hzeta=0.0;

  if (r==0.0)
    {
      hzeta=1.0;
    }
  else
    {
      hatx=x/r;
      s=boostublas::prec_inner_prod(k,hatx);
      //In case of roundoff errors
      if (s<-1.0)
	{
	  s=-1.0;
	}
      if (s>1.0)
	{
	  s=1.0;
	}
      //
      if (s!=1.0) //s not equal to 1 -> compute integral
  	{
	  gsl_integration_qng(&F,s,1.0,epsabs,epsrel,&hzeta,&error,&neval);
  	}
      hzeta+=gsl_sf_exprel(-kappar);
    }
  hzeta*=-kappa/(4.0*Pi);
}
