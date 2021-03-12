#include "pcc.h"
#include "dev_to_val.h"

/*
  Computes the values of the densities at quadrature points,
  given their spherical harmonics series:
  f(q(z)) = sum f_n^m Y_n^m(z)
  where q is the map from the sphere to the surface.
*/

void block_dev_to_val(MatComplex &valdens,MatComplex &devdens,VecConductivity &sigma)
{
  std::cout << "Computing densities values..." << std::endl;

  unsigned nsig=sigma.size();

  std::cout << nsig << " blocks" << std::endl;

  unsigned ndponetwo=devdens.size2();

  unsigned istart,iend;
  unsigned ni;

  istart=0;
  for (unsigned i=0;i<nsig;i++)
    {
      ni=pow(sigma(i).GetMaxDegree()+1,2);
      iend=istart+ni;
      std::cout << "Block " << i << std::endl;
      MatComplex devdensi = MatRangeComplex(devdens,boostublas::range(istart,iend),boostublas::range(0,ndponetwo));
      MatComplex valdensi (2*ni,ndponetwo);
      dev_to_val(valdensi,devdensi,sigma(i));
      MatRangeComplex(valdens,boostublas::range(2*istart,2*iend),boostublas::range(0,ndponetwo)) = valdensi;
      istart=iend;
    }
}
