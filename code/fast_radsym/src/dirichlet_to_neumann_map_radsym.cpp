#include "radsym.h"

/*
  Computes the Dirichlet-to-Neumann map of the equation
  -div sigma (nabla sigma) = 0 for a radially symmetric conductivity sigma
  using the eigenvalues for the Dirichlet-to-Neumann map, easily computable
  (see function dnmap_eigenvalues and dnmap_eigenvalues_pcc).
  The computations are similar to those for the Dirichlet-to-Neumann map for
  the equation -Delta u = 0, except that the eigenvalues "p" are replaced
  by the eigenvalues v(p).

  A function on the sphere which value is 1 at a quadrature point x and zero at the others is approximated by:
  f \approx wtt(x) \sum_{n=0}^{N-1} \sum_{m=-n}^{n} Y_n^{-m}(x) Y_n^m.
  Where wtt is the weight at x, it is given by wtt=\pi/N wt(theta_x) where wt are the weights given by the function gauss_legendre.

  The quadrature points are given by:
  x=(sqrt(1-t^2)*cos(phi),sqrt(1-t^2)*sin(phi),t).
  Where t are the N values given by the function gauss_legendre and phi=k*pi/N with 0<=k<=2*N-1. 

  The quadrature points are organised by first increasing the index of t and then the one of phi. 
*/

void dirichlet_to_neumann_map_radsym(MatDouble &dnmap,MatDouble &lf,MatDouble &cp,VecDouble &wt,VecDouble &eigs)

{
  std::cout << "Computing the Dirichlet-to-Neumann map for a radially symmetric conductivity... " << std::endl;

  unsigned i,j;
  unsigned n=wt.size();
  VecDouble wtt=Pi/((double)n)*wt;

  dnmap = ZeroMatDouble (dnmap.size1(),dnmap.size2());

#pragma omp parallel default(none) shared(dnmap,n,lf,cp,wtt,eigs) private(i,j)
  {
#pragma omp for collapse(4)
    for (unsigned ip=0;ip<2*n;ip++) //Loop on values of phi for the rows
      {
	for (unsigned it=0;it<n;it++) //Loop on values of theta for the rows
	  {
	    for (unsigned jp=0;jp<2*n;jp++) //Loop on values of phi for the columns
	      {
		for (unsigned jt=0;jt<n;jt++) //Loop on values of theta for the columns
		  {
		    i=n*ip+it; //Row of dnmap
		    j=n*jp+jt; //Column of dnmap

		    //Computing the sum
		    for (int nn=1;nn<n;nn++) //Loop on degrees of spherical harmonics (0 for nn=0)
		      {
			for (int mm=-nn;mm<=nn;mm++) //Loop on orders of spherical harmonics
			  {
			    dnmap(i,j)+=eigs(nn)*wtt(jt)*lf(jt,(nn*(nn+1))/2+abs(mm))*lf(it,(nn*(nn+1))/2+abs(mm))*cp(abs(int(ip-jp)),abs(mm));
			  }
		      }
		  }
	      }
	  }
      }
  }
}
