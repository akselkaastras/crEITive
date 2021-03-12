#include "eit.h"

/*
  Computes the Single Layer Operator Matrix for the quadrature points are given by:
  x=(sqrt(1-t^2)*cos(phi),sqrt(1-t^2)*sin(phi),t).
  Where t are the N values given by the function gauss_legendre and phi=k*pi/N with 0<=k<=2*N-1.

  The quadrature points are organised by first increasing the index of t and then the one of phi.

  SLO is a 2N^2x2N^2 matrix.
  It is computed according to the quadrature rule (Wienert, Kress)
  1/(4*Pi) \int_{S^2} f(y)/|x-y| ds(y) ~ 1/(4*Pi) \sum_{k=1}^N \sum_{l=0}^{2N-1} wt(x_{kl}) f(x_{kl}) \sum_{m=0}^{N-1} P_m(x_{kl}.x)
  Where wt is the weight at x_{kl}, it is given by wt=\pi/N wtl(theta_{_{kl}}) where wtl are the weights given by the function gauss_legendre and where P_m are the Legendre polynomials.
*/

void slo_matrix(MatDoubleColMaj &SLO,MatDouble &x,VecDouble &wtrep)
{
  std::cout << "Computing the Single Layer Operator Matrix... " << std::endl;

  unsigned tn2=x.size1();
  unsigned n=intsqrt(tn2/2);

#pragma omp parallel default(none) shared(SLO,x,wtrep,tn2,n)
  {
#pragma omp for collapse(2)
    for (unsigned j=0;j<tn2;j++) //Column of SLO
      {
	for (unsigned i=0;i<tn2;i++) //Row of SLO
	  {
	    double xidotxj=boostublas::prec_inner_prod(boostublas::row(x,i),boostublas::row(x,j));
	    double plm2=1.0;
	    double plm1=xidotxj;
	    SLO(i,j)=plm2+plm1;
	    for (unsigned l=2;l<n;l++) //Sum of P_l(xi.xj)
	      {
		double pl=boostmath::legendre_next(l-1,xidotxj,plm1,plm2);
		SLO(i,j)+=pl;
		plm2=plm1;
		plm1=pl;
	      }
	    SLO(i,j)*=wtrep(j);
	  }
      }
  }
  //Divide by 4*Pi
  SLO/=4.0*Pi;
}
