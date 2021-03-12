#include "pcc.h"

/*
  Computes the matrix of the operator Kbar':
  (Kbar'\varphi)(z) = \int_{surface_b} \frac{\partial\bar\Phi(z,z')}{\partial \nu(z)}
  for z \in surface_a. surface_a and surface_b can be the same.
  \bar\Phi(z,z') = 1 / ( (4*Pi) |z'| | z - z'/|z'|^2 | )
  See 3.6 of "Inverse Acoustic and Electromagnetic Scattering Theory" by
  David Colton and Rainer Kress.
*/

void kbarp_matrix(MatComplex &Kbarp,Conductivity &sigmaa,Conductivity &sigmab)
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
  double normxni2,normxnj2;
  MatDouble F (2*naponetwo,2*nbponetwo);

  //compute kernel F
  for (unsigned ip=0;ip<2*(na+1);ip++) 
    {
      for (unsigned it=0;it<na+1;it++) 
	{
	  ni=it+(na+1)*ip;
	  xni=boostublas::row(xa,ni);
	  normxni2=pow(xni(0),2)+pow(xni(1),2)+pow(xni(2),2);
	  for (unsigned jp=0;jp<2*(nb+1);jp++) 
	    {
	      for (unsigned jt=0;jt<nb+1;jt++) 
		{
		  nj=jt+(nb+1)*jp;
		  xnj=boostublas::row(xb,nj);
		  normxnj2=pow(xnj(0),2)+pow(xnj(1),2)+pow(xnj(2),2);
		  diffx=normxnj2*xni-xnj;
		  F(ni,nj)=-alphaa(it)*alphab(jt)*jyb(nj)*(boostublas::prec_inner_prod(diffx,boostublas::row(nua,ni)))/pow(sqrt(1.0-2.0*boostublas::prec_inner_prod(xni,xnj)+normxni2*normxnj2),3);
		}
	    }
	}
    }

  Kbarp = ZeroMatComplex (Kbarp.size1(),Kbarp.size2());

  //compute Kbarp
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
			  Kbarp(kl,klprime)+=F(ni,nj)*lfa(it,(k*(k+1))/2+abs(l))*lfb(jt,(kprime*(kprime+1))/2+abs(lprime))*epa(ip,na-l)*epb(jp,nb+lprime);
			}
		    }
		}
	    }
	}
    }
  Kbarp/=4.0*Pi;
}
