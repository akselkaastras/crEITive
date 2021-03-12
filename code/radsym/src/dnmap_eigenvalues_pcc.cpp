#include "radsym.h"

/*
  Computes nd+1 eigenvalues of the Dirichlet-to-Neumann map for a
  piecewise constant conductivity with values
  sigma(k) in (r(k-1),r(k)), k=0,...,n.
  sigma is a vector of size n and r a vector of size n-1.
  By convention, in the previous notation, r(-1)=0 and r(n)=1.

  The eigenvalues are given by:
  lambda(0) = 0
  lambda(p) = p-(2p+1)/(1+c(p,n-2)) for p>0
  where for k=0..n-2
  C(p,k) = w(p,k)*(b(p)*sigma(k+1)*a(p,k)+sigma(k))/(sigma(k+1)*a(p,k)-sigma(k))
  with
  a(p,0)=1
  a(p,k)=(C(p,k-1)+w(p,k))/(C(p,k-1)-b(p)*w(p,k)) for k>0
  w(p,k)=r(k)^(-(2*p+1)) for k>=0
  b(p)=(p+1)/p
*/

void dnmap_eigenvalues_pcc(VecDouble &eigenvalues,VecDouble &sigma,VecDouble &r)
{
  unsigned n=r.size();
  unsigned nd=eigenvalues.size();

  if (n==0) //Constant conductivity on [0,1]
    {
      eigenvalues(0)=0.0;
      for (unsigned p=1;p<nd;p++)
	{
	  eigenvalues(p)=(double)p*sigma(0);
	}
    }
  else
    {
      VecDouble b (nd-1);
      double a,c,w;

      for (unsigned p=1;p<nd;p++)
	{
	  b(p-1)=((double)p+1.0)/(double)p;
	}

      eigenvalues(0)=0;

      for (unsigned p=1;p<nd;p++)
	{
	  w=1.0/pow(r(0),2*p+1);
	  c=w*(b(p-1)*sigma(1)+sigma(0))/(sigma(1)-sigma(0));
	  for (unsigned k=1;k<n-1;k++)
	    {
	      w=1.0/pow(r(k),2*p+1);
	      a=(c+w)/(c-b(p-1)*w);
	      c=w*(b(p-1)*sigma(k+1)*a+sigma(k))/(sigma(k+1)*a-sigma(k));
	    }
	  eigenvalues(p)=(double)p-(2.0*(double)p+1.0)/(1.0+c);
	}
    }
}
