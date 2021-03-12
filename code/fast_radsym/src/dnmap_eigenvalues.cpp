#include "radsym.h"
#include "conductivity_on_radius.h"
#include "dnmap_eigenvalues_pcc.h"

/*
  Computes the eigenvalues of the Dirichlet-to-Neumann map using
  piecewise constant conductivity bounds sigmadown <= sigma <= sigmaup.
  Eigenvalues of Dirichlet-to-Neumann for piecewise constant conductivities
  is easily computable (see dnmap_eigenvalues_pcc) and the eigenvalues lambda(p)
  of the Dirichlet-to-Neumann for sigma verify
  lambdadown(p) <= lambda (p) <= lambdaup(p)
  where lambdadown(p), lambdaup(p) are the eigenvalues of the
  Dirichlet-to-Neumann map for respectively lambdadown and lambdaup.
*/

void dnmap_eigenvalues(VecDouble &eigenvalues,VecDouble &eigenvalueserror,MatDouble &sigma,double &stepsize)
{
  std::cout << "Computing the Dirichlet-to-Neumann map eigenvalues... " << std::endl;

  unsigned nd=eigenvalues.size()-1;
  VecDouble eigenvaluesup (nd+1);
  VecDouble eigenvaluesdown (nd+1);

  unsigned n=ceil(1.0/stepsize)+1;
  VecDouble valsigma (n);
  VecDouble r (n);
  VecDouble sigmaup (n-1);
  VecDouble rup (n-2);
  VecDouble sigmadown (n-1);
  VecDouble rdown (n-2);
  unsigned nup,ndown;

  //Conductivity at discretization points
  conductivity_on_radius(r,valsigma,sigma);

  //sigmadown <= valsigma <= sigmaup
  for (unsigned m=0;m<n-1;m++)
    {
      if (valsigma(m+1)>valsigma(m))
	{
	  sigmaup(m)=valsigma(m+1);
	  sigmadown(m)=valsigma(m);
	}
      else
	{
	  sigmaup(m)=valsigma(m);
	  sigmadown(m)=valsigma(m+1);
	}
    }

  //Avoid adjacent equal values
  nup=0;
  ndown=0;
  for (unsigned m=1;m<n-1;m++)
    {
      if (fabs(sigmaup(m)/sigmaup(nup)-1.0)>eps)
	{
	  nup+=1;
	  sigmaup(nup)=sigmaup(m);
	  rup(nup-1)=r(m);
	}
      if (fabs(sigmadown(m)/sigmadown(ndown)-1.0)>eps)
	{
	  ndown+=1;
	  sigmadown(ndown)=sigmadown(m);
	  rdown(ndown-1)=r(m);
	}
    }
  sigmaup.resize(nup+1);
  rup.resize(nup);
  sigmadown.resize(ndown+1);
  rdown.resize(ndown);

  dnmap_eigenvalues_pcc(eigenvaluesup,sigmaup,rup);
  dnmap_eigenvalues_pcc(eigenvaluesdown,sigmadown,rdown);

  eigenvalues=(eigenvaluesup+eigenvaluesdown)/2.0;
  eigenvalueserror(0)=0.0;
  for (unsigned p=1;p<=nd;p++)
    {
      eigenvalueserror(p)=(eigenvaluesup(p)-eigenvaluesdown(p))/(fabs(eigenvalues(p)-(double)p));
    }
}

