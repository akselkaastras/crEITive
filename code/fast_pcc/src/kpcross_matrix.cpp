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
			    ni=it+(na+1)*ip;
			    nj=jt+(nb+1)*jp;
			    xni=boostublas::row(xa,ni);
			    xnj=boostublas::row(xb,nj);
			    diffx=xni-xnj;
			    diffx/=pow(boostublas::norm_2(diffx),3);
			    F(ni,nj)=-alphaa(it)*alphab(jt)*jyb(nj)*(boostublas::prec_inner_prod(diffx,boostublas::row(nua,ni)));
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
			    ni=it+(na+1)*ip; // quadrature points on surface a
			    nj=jt+(nb+1)*jp; // quadrature points on surface b
			    xni=boostublas::row(xa,ni);
			    xnj=boostublas::row(xb,nj);
			    xnj2=boostublas::row(xb2,nj);
			    diffx=xni-xnj;
			    diffx2=xni-xnj2;
			    diffx/=pow(boostublas::norm_2(diffx),3);
			    diffx2/=pow(boostublas::norm_2(diffx2),3);
			    F(ni,nj)=-alphaa(it)*alphab(jt)*jyb(nj)*(boostublas::prec_inner_prod(diffx,boostublas::row(nua,ni)));
			    Fb(ni,nj)=-alphaa(it)*alphab(jt)*jyb2(nj)*(boostublas::prec_inner_prod(diffx2,boostublas::row(nua,ni)));
			  }
		      }
		  }
	      }
	  }
  }
  else if (widthb > 0.01 && widtha >0.01)
  {


  	#pragma omp parallel default(none) shared(F,Fa,Fb,Fab,xa,xb,xa2,xb2,alphaa,alphab,jyb,jyb2,nua,nua2,na,nb) private(ni,nj)
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
			    ni=it+(na+1)*ip; // quadrature points on surface a
			    nj=jt+(nb+1)*jp; // quadrature points on surface b
			    xni=boostublas::row(xa,ni);
			    xnj=boostublas::row(xb,nj);
			    xni2=boostublas::row(xa2,ni);
			    xnj2=boostublas::row(xb2,nj);
			    diffx=xni-xnj;
			    diffx2=xni2-xnj;
			    diffx3=xni-xnj2;
			    diffx4=xni2-xnj2;
			    diffx/=pow(boostublas::norm_2(diffx),3);
			    diffx2/=pow(boostublas::norm_2(diffx2),3);
			    diffx3/=pow(boostublas::norm_2(diffx3),3);
			    diffx4/=pow(boostublas::norm_2(diffx4),3);
			    F(ni,nj)=-alphaa(it)*alphab(jt)*jyb(nj)*(boostublas::prec_inner_prod(diffx,boostublas::row(nua,ni)));
			    Fa(ni,nj)=-alphaa(it)*alphab(jt)*jyb(nj)*(boostublas::prec_inner_prod(diffx2,boostublas::row(nua2,ni)));
			    Fb(ni,nj)=-alphaa(it)*alphab(jt)*jyb2(nj)*(boostublas::prec_inner_prod(diffx3,boostublas::row(nua,ni)));
			    Fab(ni,nj)=-alphaa(it)*alphab(jt)*jyb2(nj)*(boostublas::prec_inner_prod(diffx4,boostublas::row(nua2,ni)));
			  }
		      }
		  }
	      }
	  }
  }
  else
  {

  	#pragma omp parallel default(none) shared(F,Fa,xa,xb,xa2,alphaa,alphab,jyb,nua,nua2,na,nb) private(ni,nj)
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
			    ni=it+(na+1)*ip; // quadrature points on surface a
			    nj=jt+(nb+1)*jp; // quadrature points on surface b
			    xni=boostublas::row(xa,ni);
			    xnj=boostublas::row(xb,nj);
			    xni2=boostublas::row(xa2,ni);
			    diffx=xni-xnj;
			    diffx2=xni2-xnj;
			    diffx/=pow(boostublas::norm_2(diffx),3);
			    diffx2/=pow(boostublas::norm_2(diffx2),3);
			    F(ni,nj)=-alphaa(it)*alphab(jt)*jyb(nj)*(boostublas::prec_inner_prod(diffx,boostublas::row(nua,ni)));
			    Fa(ni,nj)=-alphaa(it)*alphab(jt)*jyb(nj)*(boostublas::prec_inner_prod(diffx2,boostublas::row(nua2,ni)));
			  }
		      }
		  }
	      }
	  }
  }


  Kpcross = ZeroMatComplex (Kpcross.size1(),Kpcross.size2());



if (widthb < 0.01 && widtha < 0.01)
  {
  	  //compute Kpcross
	#pragma omp parallel default(none) shared(Kpcross,F,lfa,epa,lfb,epb,na,nb,naponetwo,nbponetwo) private(k,l,kprime,lprime,ni,nj)
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
				    Kpcross(kl,klprime)+=F(ni,nj)*lfa(it,(k*(k+1))/2+abs(l))*lfb(jt,(kprime*(kprime+1))/2+abs(lprime))*epa(ip,na-l)*epb(jp,nb+lprime);
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
  	 //compute Kpcross
	#pragma omp parallel default(none) shared(Kpcross,F,Fb,lfa,epa,lfb,epb,na,nb,naponetwo,nbponetwo) private(k,l,kprime,lprime,ni,nj)
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
				    Kpcross(kl,klprime)+=F(ni,nj)*lfa(it,(k*(k+1))/2+abs(l))*lfb(jt,(kprime*(kprime+1))/2+abs(lprime))*epa(ip,na-l)*epb(jp,nb+lprime);
				    Kpcross(kl,klprime)+=Fb(ni,nj)*lfa(it,(k*(k+1))/2+abs(l))*lfb(jt,(kprime*(kprime+1))/2+abs(lprime))*epa(ip,na-l)*epb(jp,nb+lprime);
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
   	 //compute Kpcross
	#pragma omp parallel default(none) shared(Kpcross,F,Fa,Fb,Fab,lfa,epa,lfb,epb,na,nb,naponetwo,nbponetwo) private(k,l,kprime,lprime,ni,nj)
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
				    Kpcross(kl,klprime)+=F(ni,nj)*lfa(it,(k*(k+1))/2+abs(l))*lfb(jt,(kprime*(kprime+1))/2+abs(lprime))*epa(ip,na-l)*epb(jp,nb+lprime);
				    Kpcross(kl,klprime)+=Fb(ni,nj)*lfa(it,(k*(k+1))/2+abs(l))*lfb(jt,(kprime*(kprime+1))/2+abs(lprime))*epa(ip,na-l)*epb(jp,nb+lprime);
				    Kpcross(kl,klprime)+=Fa(ni,nj)*lfa(it,(k*(k+1))/2+abs(l))*lfb(jt,(kprime*(kprime+1))/2+abs(lprime))*epa(ip,na-l)*epb(jp,nb+lprime);
				    Kpcross(kl,klprime)+=Fab(ni,nj)*lfa(it,(k*(k+1))/2+abs(l))*lfb(jt,(kprime*(kprime+1))/2+abs(lprime))*epa(ip,na-l)*epb(jp,nb+lprime);
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
  	 //compute Kpcross
	#pragma omp parallel default(none) shared(Kpcross,F,Fa,lfa,epa,lfb,epb,na,nb,naponetwo,nbponetwo) private(k,l,kprime,lprime,ni,nj)
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
				    Kpcross(kl,klprime)+=F(ni,nj)*lfa(it,(k*(k+1))/2+abs(l))*lfb(jt,(kprime*(kprime+1))/2+abs(lprime))*epa(ip,na-l)*epb(jp,nb+lprime);
				    Kpcross(kl,klprime)+=Fa(ni,nj)*lfa(it,(k*(k+1))/2+abs(l))*lfb(jt,(kprime*(kprime+1))/2+abs(lprime))*epa(ip,na-l)*epb(jp,nb+lprime);
				  }
			      }
			  }
		      }
		  }
	      }
	  }
  }


  Kpcross/=4.0*Pi;


}
