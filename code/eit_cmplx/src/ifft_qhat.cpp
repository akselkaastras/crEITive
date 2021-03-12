#include "eit.h"

/*
  Computes the Inverse Fourier Transform of qhat. 
*/

void ifft_qhat(MatDouble &xi,Array3Complex &qhat)
{
  std::cout << "Computing the Inverse Fourier Transform of qhat..." << std::endl;

  unsigned ng=qhat.shape()[0];
  unsigned ng2=square(ng);
  unsigned ng3=ng*ng2;
  unsigned i,j,k;
  double nmoon=((double)ng-1.0)/(double)ng;
  double ximax=-xi(0,0);

  // Phase
  VecDouble anglephase (ng);
  std::copy(CountItDouble0,CountItDouble0+ng,anglephase.begin());
  anglephase*=-nmoon*Pi;
  VecComplex phase = complexfunvec(VecOp(anglephase,[](double c) -> double { return cos(c); }),VecOp(anglephase,[](double c) -> double { return sin(c); }));

  // Dephasing
  for (unsigned p=0;p<ng3;p++)
    {
      i=p/ng2;
      j=(p/ng)%ng; //Same as (p%ng2)/ng
      k=p%ng; //Same as (p%ng2)%ng
      qhat[i][j][k]*=phase(i)*phase(j)*phase(k);
    }

  // Computing ifft
  fftw_plan plan;
  plan=fftw_plan_dft_3d(ng,ng,ng,reinterpret_cast<fftw_complex *>(&qhat[0][0][0]),reinterpret_cast<fftw_complex *>(&qhat[0][0][0]),FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  // Dephasing and multiplying by coefficient
  dcomplex coeff=pow(nmoon,3)/8.0*dcomplex(cos(3.0*ximax),sin(3.0*ximax));
  for (unsigned p=0;p<ng3;p++)
    {
      i=p/ng2;
      j=(p/ng)%ng; //Same as (p%ng2)/ng
      k=p%ng; //Same as (p%ng2)%ng
      qhat[i][j][k]*=coeff*phase(i)*phase(j)*phase(k);
    }
}
