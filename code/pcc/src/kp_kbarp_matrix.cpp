#include "pcc.h"
#include "kp_matrix.h"
#include "kpcross_matrix.h"
#include "kbarp_matrix.h"

/*
  Computes the matrix of the operator I-2 mu(K'-Kbar'):
  (K'\varphi)(z) = \int_{surface} \frac{\partial\Phi(z,z')}{\partial \nu(z)}
  (Kbar'\varphi)(z) = \int_{surface} \frac{\partial\bar\Phi(z,z')}{\partial \nu(z)}
  \bar\Phi(z,z') = 1 / ( (4*Pi) |z'| | z - z'/|z'|^2 | )
  See 3.6 of "Inverse Acoustic and Electromagnetic Scattering Theory" by
  David Colton and Rainer Kress.
*/

void kp_kbarp_matrix(MatComplex &KpKbarp,VecConductivity &sigma)
{
  std::cout << "Building matrix..." << std::endl;

  unsigned nsig=sigma.size();

  std::cout << nsig << " x " << nsig << " blocks" << std::endl;

  double twomu;
  unsigned istart,jstart,iend,jend;
  unsigned ni,nj;

  istart=0;
  for (unsigned i=0;i<nsig;i++)
    {
      ni=pow(sigma(i).GetMaxDegree()+1,2);
      iend=istart+ni;
      twomu=sigma(i).GetTwoMu();
      jstart=0;
      for (unsigned j=0;j<nsig;j++)
	{
	  std::cout << "Block (" << i << "," << j << ")..." << std::endl;
	  nj=pow(sigma(j).GetMaxDegree()+1,2);
	  jend=jstart+nj;
	  MatComplex Kp (ni,nj);
	  MatComplex Kbarp (ni,nj);
	  if (i==j)
	    {
	      kp_matrix(Kp,sigma(i));
	    }
	  else
	    {
	      kpcross_matrix(Kp,sigma(i),sigma(j));
	    }
	  kbarp_matrix(Kbarp,sigma(i),sigma(j));
	  MatRangeComplex(KpKbarp,boostublas::range(istart,iend),boostublas::range(jstart,jend)) = -twomu*(Kp-Kbarp);
	  jstart=jend;
	}
      istart=iend;
    }

  // Add identity matrix
  for (unsigned ij=0;ij<KpKbarp.size2();ij++)
    {
      KpKbarp(ij,ij)+=1.0;
    }
}
