#include "pcc.h"

/*
  Computes the matrix of the operator K':
  (K'\varphi)(z) = \int_{surface_b} \frac{\partial\Phi(z,z')}{\partial \nu(z)}
  For z in surface_a, where surface_a and surface_b are distinct.
  See 3.6 of "Inverse Acoustic and Electromagnetic Scattering Theory" by
  David Colton and Rainer Kress.
*/

void kpcross_matrix(MatComplex &Kpcross,Conductivity &sigmaa,Conductivity &sigmab)
{
  unsigned na=sigmaa.GetMaxDegree();
  VecDouble alphaa=sigmaa.GetAlpha();
  MatDouble xa=sigmaa.GetSurfacePoints();
  MatDouble nua=sigmaa.GetOutwardNormal();
  MatDouble lfa=sigmaa.GetAssociatedLegendreFunctions();
  MatComplex epa=sigmaa.GetComplexExponentials();

  unsigned nb=sigmab.GetMaxDegree();
  VecDouble alphab=sigmab.GetAlpha();
  MatDouble xb=sigmab.GetSurfacePoints();
  VecDouble jyb=sigmab.GetJacobian();
  MatDouble lfb=sigmab.GetAssociatedLegendreFunctions();
  MatComplex epb=sigmab.GetComplexExponentials();

  unsigned naponetwo=pow(na+1,2);
  unsigned nbponetwo=pow(nb+1,2);
  unsigned ni,nj,k,kprime;
  int l,lprime;

  VecDouble xni (3);
  VecDouble xnj (3);
  VecDouble diffx (3);
  MatDouble F (2*naponetwo,2*nbponetwo);

  //compute kernel F
  for (unsigned ip=0;ip<2*(na+1);ip++) 
    {
      for (unsigned it=0;it<na+1;it++) 
	{
	  ni=it+(na+1)*ip;
	  xni=boostublas::row(xa,ni);
	  for (unsigned jp=0;jp<2*(nb+1);jp++) 
	    {
	      for (unsigned jt=0;jt<nb+1;jt++) 
		{
		  nj=jt+(nb+1)*jp;
		  xnj=boostublas::row(xb,nj);
		  diffx=xni-xnj;
		  diffx/=pow(boostublas::norm_2(diffx),3);
		  F(ni,nj)=-alphaa(it)*alphab(jt)*jyb(nj)*(boostublas::prec_inner_prod(diffx,boostublas::row(nua,ni)));
		}
	    }
	}
    }

  Kpcross = ZeroMatComplex (Kpcross.size1(),Kpcross.size2());

  //compute Kpcross
  for (unsigned kl=0;kl<naponetwo;kl++) 
    {
      k=intsqrt(kl);//degree
      l=kl-k*k-k;//order
      for (unsigned klprime=0;klprime<nbponetwo;klprime++) 
	{
	  kprime=intsqrt(klprime);//degree
	  lprime=klprime-kprime*kprime-kprime;//order

	  //make the sum
	  for (unsigned ip=0;ip<2*(na+1);ip++) 
	    {
	      for (unsigned it=0;it<na+1;it++) 
		{
		  ni=it+(na+1)*ip;
		  for (unsigned jp=0;jp<2*(nb+1);jp++) 
		    {
		      for (unsigned jt=0;jt<nb+1;jt++) 
			{
			  nj=jt+(nb+1)*jp;
			  Kpcross(kl,klprime)+=F(ni,nj)*lfa(it,(k*(k+1))/2+abs(l))*lfb(jt,(kprime*(kprime+1))/2+abs(lprime))*epa(ip,na-l)*epb(jp,nb+lprime);
			}
		    }
		}
	    }
	}
    }
  Kpcross/=4.0*Pi;
}
