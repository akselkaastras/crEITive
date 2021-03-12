#include "eit.h"

/*
  Computes the Dirichlet-to-Neumann map of the equation -Delta u=0.

  A function on the sphere which value is 1 at a quadrature point x and zero at the others is approximated by:
  f \approx wtt(x) \sum_{n=0}^{N-1} \sum_{m=-n}^{n} Y_n^{-m}(x) Y_n^m.
  Where wtt is the weight at x, it is given by wtt=\pi/N wt(theta_x) where wt are the weights given by the function gauss_legendre.

  The quadrature points are given by:
  x=(sqrt(1-t^2)*cos(phi),sqrt(1-t^2)*sin(phi),t).
  Where t are the N values given by the function gauss_legendre and phi=k*pi/N with 0<=k<=2*N-1. 

  The quadrature points are organised by first increasing the index of t and then the one of phi. 
*/

void dirichlet_to_neumann_map_0(MatDoubleColMaj &dnmap,MatDouble &lf,MatDouble &cp,VecDouble &wt)
{
  std::cout << "Computing the Dirichlet-to-Neumann map for a conductivity equal to 1... " << std::endl;

  unsigned n=wt.size();
  unsigned tn2=2*square(n);
  unsigned vn2=(n*(n+1))/2;
  unsigned nn,mm;
  unsigned in,jn;
  unsigned ip,jp;
  dnmap = ZeroMatDoubleColMaj (dnmap.size1(),dnmap.size2());

  for (unsigned j=0;j<tn2;j++) //Column of dnmap
    {
      jn=j%n;
      jp=j/n;
      for (unsigned i=0;i<tn2;i++) //Row of dnmap
	{
	  in=i%n;
	  ip=i/n;
	  //Computing the sum
	  for (unsigned nm=1;nm<vn2;nm++) //Loop on degrees and orders of spherical harmonics
	    {
	      nn=(intsqrt(1+8*nm)-1)/2; //Degree of spherical harmonic
	      mm=nm-(nn*(nn+1))/2; //Order of spherical harmonic
	      dnmap(i,j)+=onetwo(mm)*(double)nn*wt(jn)*lf(jn,nm)*lf(in,nm)*cp(abs((int)ip-(int)jp),mm);
	    }
	}
    }
}

