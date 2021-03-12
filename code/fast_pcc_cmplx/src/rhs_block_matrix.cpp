#include "pcc.h"
#include "rhs_matrix.h"

/*
  Computes the matrix of the right hand side:
  2 mu \int_S (\partial_\nu R_{n'}^{m'}) (q(z)) Y_n^{-m} (z) ds(z)
  Where R_{n'}^{m'} are the solid harmonics, Y_n^{-m} the spherical harmonics,
  q is the map from the sphere S to the surface.
*/

void rhs_block_matrix(MatComplex &rhs,VecConductivity &sigma,unsigned &nd)
{
  std::cout << "Building right hand side..." << std::endl;

  unsigned nsig=sigma.size();

  std::cout << nsig << " blocks" << std::endl;

  unsigned ndponetwo=(nd+1)*(nd+1);

  dcomplex twomu;
  unsigned istart,iend;
  unsigned ni;

  istart=0;
  for (unsigned i=0;i<nsig;i++)
    {
      ni=pow(sigma(i).GetMaxDegree()+1,2);
      iend=istart+ni;
      twomu=sigma(i).GetTwoMu();
      std::cout << "Block " << i << std::endl;
      MatComplex rhsi (ni,ndponetwo);
      rhs_matrix(rhsi,sigma(i),nd);
      MatRangeComplex(rhs,boostublas::range(istart,iend),boostublas::range(0,ndponetwo)) = twomu*rhsi;
      istart=iend;
    }
}
