#include "pcc.h"

/*
  Computes the Dirichlet-to-Neumann map of the equation -div(\sigma grad u)=0.

  A function on the sphere which value is 1 at a quadrature point x and zero at the others is approximated by:
  f \approx alpha(x) \sum_{n=0}^{N} \sum_{m=-n}^{n} Y_n^{-m}(x) Y_n^m.
  Where wtt is the weight at x, it is given by alpha=\pi/(N+1) w(theta_x) where w are the weights given by the function gauss_legendre.

  The quadrature points are given by:
  x=(sqrt(1-t^2)*cos(phi),sqrt(1-t^2)*sin(phi),t).
  Where t are the N values given by the function gauss_legendre and phi=k*pi/(N+1) with 0<=k<=2*N+1. 

  The quadrature points are organised by first increasing the index of t and then the one of phi. 
*/

void dirichlet_to_neumann_map(MatComplex &DNmap,MatComplex &ND,MatDouble &lf,MatComplex &ep,VecDouble &alpha)
{
  std::cout << "Computing the Dirichlet-to-Neumann map... " << std::endl;

  unsigned i,j;
  unsigned k;
  int l;
  unsigned npone=alpha.size();

  DNmap = ZeroMatComplex (DNmap.size1(),DNmap.size2());

  for (unsigned ip=0;ip<2*npone;ip++) //Rows, phi
    {
      for (unsigned it=0;it<npone;it++) //Rows, theta
	{
	  i=it+npone*ip; //Row of DNmap
	  for (unsigned jp=0;jp<2*npone;jp++) //Columns, phi
	    {
	      for (unsigned jt=0;jt<npone;jt++) //Columns, theta
		{
		  j=jt+npone*jp; //Column of DNmap

		  //Computing the sum
		  for (unsigned kl=0;kl<npone*npone;kl++)//degree and order of spherical harmonics
		    {
		      k=intsqrt(kl);//degree
		      l=kl-k*k-k;//order
		      DNmap(i,j)+=alpha(jt)*lf(jt,(k*(k+1))/2+abs(l))*ep(jp,npone-1-l)*(ND(i,k*(k+1)+l)+(double)k*lf(it,(k*(k+1))/2+abs(l))*ep(ip,npone-1+l));
		    }
		}
	    }
	}
    }
}
