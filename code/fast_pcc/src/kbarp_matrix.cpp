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
  double widtha = sigmaa.GetWidth();

  unsigned nb=sigmab.GetMaxDegree();
  VecDouble alphab=sigmab.GetAlpha();
  MatDouble xb=sigmab.GetSurfacePoints();
  VecDouble jyb=sigmab.GetJacobian();
  MatDouble lfb=sigmab.GetAssociatedLegendreFunctions();
  MatComplex epb=sigmab.GetComplexExponentials();
  double widthb = sigmab.GetWidth();

  unsigned naponetwo=pow(na+1,2);
  unsigned nbponetwo=pow(nb+1,2);
  unsigned ni,nj,k,kprime;
  int l,lprime;

  MatDouble F (2*naponetwo,2*nbponetwo);
  MatDouble xb2=sigmab.GetSurfacePoints2();
  VecDouble jyb2=sigmab.GetJacobian2();
  MatDouble Fb (2*naponetwo,2*nbponetwo);
  MatDouble xa2=sigmaa.GetSurfacePoints2();
  MatDouble nua2=sigmaa.GetOutwardNormal2();
  MatDouble Fa (2*naponetwo,2*nbponetwo);
  MatDouble Fab (2*naponetwo,2*nbponetwo);

  //compute kernel F
  if (widthb < 0.01 && widtha < 0.01)
  {
	#pragma omp parallel default(none) shared(F,xa,xb,alphaa,alphab,jyb,nua,na,nb) private(ni,nj)
	  {
	#pragma omp for collapse(4)
	    for (unsigned ip=0;ip<2*(na+1);ip++) 
	      {
		for (unsigned it=0;it<na+1;it++) 
		  {
		    for (unsigned jp=0;jp<2*(nb+1);jp++) 
		      {
			for (unsigned jt=0;jt<nb+1;jt++) 
			  {
			    VecDouble xni (3);
			    VecDouble xnj (3);
			    VecDouble diffx (3);
			    double normxni2,normxnj2;
			    ni=it+(na+1)*ip;
			    nj=jt+(nb+1)*jp;
			    xni=boostublas::row(xa,ni);
			    xnj=boostublas::row(xb,nj);
			    normxni2=pow(xni(0),2)+pow(xni(1),2)+pow(xni(2),2);
			    normxnj2=pow(xnj(0),2)+pow(xnj(1),2)+pow(xnj(2),2);
			    diffx=normxnj2*xni-xnj;
			    F(ni,nj)=-alphaa(it)*alphab(jt)*jyb(nj)*(boostublas::prec_inner_prod(diffx,boostublas::row(nua,ni)))/pow(sqrt(1.0-2.0*boostublas::prec_inner_prod(xni,xnj)+normxni2*normxnj2),3);
			  }
		      }
		  }
	      }
	  }
}
else if (widthb >0.01 && widtha <0.01)
{

  	#pragma omp parallel default(none) shared(F,Fb,xa,xb,xb2,alphaa,alphab,jyb,jyb2,nua,na,nb) private(ni,nj)
	  {
	#pragma omp for collapse(4)
	    for (unsigned ip=0;ip<2*(na+1);ip++) 
	      {
		for (unsigned it=0;it<na+1;it++) 
		  {
		    for (unsigned jp=0;jp<2*(nb+1);jp++) 
		      {
			for (unsigned jt=0;jt<nb+1;jt++) 
			  {
			    VecDouble xni (3);
			    VecDouble xnj (3);
			    VecDouble xnj2 (3);
			    VecDouble diffx (3);
			    VecDouble diffx2 (3);
			    double normxni2,normxnj2,normxnj22;
			    ni=it+(na+1)*ip;
			    nj=jt+(nb+1)*jp;
			    xni=boostublas::row(xa,ni);
			    xnj=boostublas::row(xb,nj);
			    xnj2=boostublas::row(xb2,nj);
			    normxni2=pow(xni(0),2)+pow(xni(1),2)+pow(xni(2),2);
			    normxnj2=pow(xnj(0),2)+pow(xnj(1),2)+pow(xnj(2),2);
			    normxnj22=pow(xnj2(0),2)+pow(xnj2(1),2)+pow(xnj2(2),2); 
			    diffx=normxnj2*xni-xnj;
			    diffx2=normxnj22*xni-xnj2;
			    F(ni,nj)=-alphaa(it)*alphab(jt)*jyb(nj)*(boostublas::prec_inner_prod(diffx,boostublas::row(nua,ni)))/pow(sqrt(1.0-2.0*boostublas::prec_inner_prod(xni,xnj)+normxni2*normxnj2),3);
			  	Fb(ni,nj)=-alphaa(it)*alphab(jt)*jyb2(nj)*(boostublas::prec_inner_prod(diffx2,boostublas::row(nua,ni)))/pow(sqrt(1.0-2.0*boostublas::prec_inner_prod(xni,xnj2)+normxni2*normxnj22),3);
			  
			  }
		      }
		  }
	      }
	  }
}
else if (widthb > 0.01 && widtha >0.01)
{
	
  	#pragma omp parallel default(none) shared(F,Fb,Fa,Fab,xa,xa2,xb,xb2,alphaa,alphab,jyb,jyb2,nua,nua2,na,nb) private(ni,nj)
	  {
	#pragma omp for collapse(4)
	    for (unsigned ip=0;ip<2*(na+1);ip++) 
	      {
		for (unsigned it=0;it<na+1;it++) 
		  {
		    for (unsigned jp=0;jp<2*(nb+1);jp++) 
		      {
			for (unsigned jt=0;jt<nb+1;jt++) 
			  {
			    VecDouble xni (3);
			    VecDouble xnj (3);
			    VecDouble xni2 (3);
			    VecDouble xnj2 (3);
			    VecDouble diffx (3);
			    VecDouble diffx2 (3);
			    VecDouble diffx3 (3);
			    VecDouble diffx4 (3);
			    double normxni2,normxnj2,normxni22,normxnj22;
			    ni=it+(na+1)*ip;
			    nj=jt+(nb+1)*jp;
			    xni=boostublas::row(xa,ni);
			    xnj=boostublas::row(xb,nj);
			    xni2=boostublas::row(xa2,ni);
			    xnj2=boostublas::row(xb2,nj);
			    normxni2=pow(xni(0),2)+pow(xni(1),2)+pow(xni(2),2);
			    normxnj2=pow(xnj(0),2)+pow(xnj(1),2)+pow(xnj(2),2);
			    normxni22=pow(xni2(0),2)+pow(xni2(1),2)+pow(xni2(2),2); 
			    normxnj22=pow(xnj2(0),2)+pow(xnj2(1),2)+pow(xnj2(2),2); 
			    diffx=normxnj2*xni-xnj;
			    diffx2=normxnj22*xni-xnj2;
			    diffx3=normxnj2*xni2-xnj;
			    diffx4=normxnj22*xni2-xnj2;
			    F(ni,nj)=-alphaa(it)*alphab(jt)*jyb(nj)*(boostublas::prec_inner_prod(diffx,boostublas::row(nua,ni)))/pow(sqrt(1.0-2.0*boostublas::prec_inner_prod(xni,xnj)+normxni2*normxnj2),3);
			  	Fa(ni,nj)=-alphaa(it)*alphab(jt)*jyb(nj)*(boostublas::prec_inner_prod(diffx3,boostublas::row(nua2,ni)))/pow(sqrt(1.0-2.0*boostublas::prec_inner_prod(xni2,xnj)+normxni22*normxnj2),3);
			  	Fb(ni,nj)=-alphaa(it)*alphab(jt)*jyb2(nj)*(boostublas::prec_inner_prod(diffx2,boostublas::row(nua,ni)))/pow(sqrt(1.0-2.0*boostublas::prec_inner_prod(xni,xnj2)+normxni2*normxnj22),3);
			  	Fab(ni,nj)=-alphaa(it)*alphab(jt)*jyb2(nj)*(boostublas::prec_inner_prod(diffx4,boostublas::row(nua2,ni)))/pow(sqrt(1.0-2.0*boostublas::prec_inner_prod(xni2,xnj2)+normxni22*normxnj22),3);
			  
			  }
		      }
		  }
	      }
	  }
}
else 
{

  	#pragma omp parallel default(none) shared(F,Fa,xa,xa2,xb,alphaa,alphab,jyb,nua,nua2,na,nb) private(ni,nj)
	  {
	#pragma omp for collapse(4)
	    for (unsigned ip=0;ip<2*(na+1);ip++) 
	      {
		for (unsigned it=0;it<na+1;it++) 
		  {
		    for (unsigned jp=0;jp<2*(nb+1);jp++) 
		      {
			for (unsigned jt=0;jt<nb+1;jt++) 
			  {
			    VecDouble xni (3);
			    VecDouble xnj (3);
			    VecDouble xni2 (3);
			    VecDouble diffx (3);
			    VecDouble diffx2 (3);
			    double normxni2,normxnj2,normxni22;
			    ni=it+(na+1)*ip;
			    nj=jt+(nb+1)*jp;
			    xni=boostublas::row(xa,ni);
			    xnj=boostublas::row(xb,nj);
			    xni2=boostublas::row(xa2,ni);
			    normxni2=pow(xni(0),2)+pow(xni(1),2)+pow(xni(2),2);
			    normxnj2=pow(xnj(0),2)+pow(xnj(1),2)+pow(xnj(2),2);
			    normxni22=pow(xni2(0),2)+pow(xni2(1),2)+pow(xni2(2),2); 
			    diffx=normxnj2*xni-xnj;
			    diffx2=normxnj2*xni2-xnj;
			    F(ni,nj)=-alphaa(it)*alphab(jt)*jyb(nj)*(boostublas::prec_inner_prod(diffx,boostublas::row(nua,ni)))/pow(sqrt(1.0-2.0*boostublas::prec_inner_prod(xni,xnj)+normxni2*normxnj2),3);
			  	Fa(ni,nj)=-alphaa(it)*alphab(jt)*jyb(nj)*(boostublas::prec_inner_prod(diffx2,boostublas::row(nua2,ni)))/pow(sqrt(1.0-2.0*boostublas::prec_inner_prod(xni2,xnj)+normxni22*normxnj2),3);
			  
			  }
		      }
		  }
	      }
	  }
}

  Kbarp = ZeroMatComplex (Kbarp.size1(),Kbarp.size2());


if (widthb < 0.01 && widtha < 0.01)
  {
  	  //compute Kpbarp
	#pragma omp parallel default(none) shared(Kbarp,F,lfa,epa,lfb,epb,na,nb,naponetwo,nbponetwo) private(k,l,kprime,lprime,ni,nj)
	  {
	#pragma omp for collapse(2) 
	    for (unsigned kl=0;kl<naponetwo;kl++) 
	      {
		for (unsigned klprime=0;klprime<nbponetwo;klprime++) 
		  {
		    k=intsqrt(kl);//degree
		    l=kl-k*k-k;//order
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
	  }
  }
  else if (widthb >0.01 && widtha <0.01)
  {
  	 //compute Kpbarp
	#pragma omp parallel default(none) shared(Kbarp,F,Fb,lfa,epa,lfb,epb,na,nb,naponetwo,nbponetwo) private(k,l,kprime,lprime,ni,nj)
	  {
	#pragma omp for collapse(2) 
	    for (unsigned kl=0;kl<naponetwo;kl++) 
	      {
		for (unsigned klprime=0;klprime<nbponetwo;klprime++) 
		  {
		    k=intsqrt(kl);//degree
		    l=kl-k*k-k;//order
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
				    Kbarp(kl,klprime)+=Fb(ni,nj)*lfa(it,(k*(k+1))/2+abs(l))*lfb(jt,(kprime*(kprime+1))/2+abs(lprime))*epa(ip,na-l)*epb(jp,nb+lprime);
				  }
			      }
			  }
		      }
		  }
	      }
	  }
  }
  else if (widthb > 0.01 && widtha >0.01)
  {
   	 //compute Kpbarp
	#pragma omp parallel default(none) shared(Kbarp,F,Fa,Fb,Fab,lfa,epa,lfb,epb,na,nb,naponetwo,nbponetwo) private(k,l,kprime,lprime,ni,nj)
	  {
	#pragma omp for collapse(2) 
	    for (unsigned kl=0;kl<naponetwo;kl++) 
	      {
		for (unsigned klprime=0;klprime<nbponetwo;klprime++) 
		  {
		    k=intsqrt(kl);//degree
		    l=kl-k*k-k;//order
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
				    Kbarp(kl,klprime)+=Fb(ni,nj)*lfa(it,(k*(k+1))/2+abs(l))*lfb(jt,(kprime*(kprime+1))/2+abs(lprime))*epa(ip,na-l)*epb(jp,nb+lprime);
				    Kbarp(kl,klprime)+=Fa(ni,nj)*lfa(it,(k*(k+1))/2+abs(l))*lfb(jt,(kprime*(kprime+1))/2+abs(lprime))*epa(ip,na-l)*epb(jp,nb+lprime);
				    Kbarp(kl,klprime)+=Fab(ni,nj)*lfa(it,(k*(k+1))/2+abs(l))*lfb(jt,(kprime*(kprime+1))/2+abs(lprime))*epa(ip,na-l)*epb(jp,nb+lprime);
				  }
			      }
			  }
		      }
		  }
	      }
	  }
  }
  else
  {
  	 //compute Kpbarp
	#pragma omp parallel default(none) shared(Kbarp,F,Fa,lfa,epa,lfb,epb,na,nb,naponetwo,nbponetwo) private(k,l,kprime,lprime,ni,nj)
	  {
	#pragma omp for collapse(2) 
	    for (unsigned kl=0;kl<naponetwo;kl++) 
	      {
		for (unsigned klprime=0;klprime<nbponetwo;klprime++) 
		  {
		    k=intsqrt(kl);//degree
		    l=kl-k*k-k;//order
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
				    Kbarp(kl,klprime)+=Fa(ni,nj)*lfa(it,(k*(k+1))/2+abs(l))*lfb(jt,(kprime*(kprime+1))/2+abs(lprime))*epa(ip,na-l)*epb(jp,nb+lprime);
				  }
			      }
			  }
		      }
		  }
	      }
	  }
  }


  Kbarp/=4.0*Pi;
}
