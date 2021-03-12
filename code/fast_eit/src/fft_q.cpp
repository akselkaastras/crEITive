//#include "moments.h"
#include "eit.h"

/*
  Computes the Fourier Transform of q and saving in same Array3Complex
*/

void fft_q(double &ximax,Array3Complex &q)
{
  std::cout << "Computing the Fourier Transform of q (for saving and comparison)..." << std::endl;

  unsigned ng=q.shape()[0];
  
  unsigned ng2=square(ng);
  unsigned ng3=ng*ng2;
  unsigned i,j,k;
  double nmoon=((double)ng-1.0)/(double)ng;

  // Phase
  VecDouble anglephase (ng);
  std::copy(CountItDouble0,CountItDouble0+ng,anglephase.begin());
  anglephase*=Pi*nmoon;
  VecComplex phase = complexfunvec(VecOp(anglephase,[](double c) -> double { return cos(c); }),VecOp(anglephase,[](double c) -> double { return sin(c); }));

  // Dephasing
  for (unsigned p=0;p<ng3;p++)
    {
      i=p/ng2;
      j=(p/ng)%ng; //Same as (p%ng2)/ng
      k=p%ng; //Same as (p%ng2)%ng
      q[i][j][k]*= phase(i)*phase(j)*phase(k);
    }

  // Computing fft
  fftw_plan plan;
  plan=fftw_plan_dft_3d(ng,ng,ng,reinterpret_cast<fftw_complex *>(&q[0][0][0]),reinterpret_cast<fftw_complex *>(&q[0][0][0]),FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  // Dephasing
  dcomplex coeff=8.0/pow(ng-1,3)*dcomplex(cos(-3.0*ximax),sin(-3.0*ximax));
  for (unsigned p=0;p<ng3;p++)
    {
      i=p/ng2;
      j=(p/ng)%ng; //Same as (p%ng2)/ng
      k=p%ng; //Same as (p%ng2)%ng
      q[i][j][k]*=coeff*phase(i)*phase(j)*phase(k);
    }
}
