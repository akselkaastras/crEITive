#include "eit.h"

/*
  Computes the conductivity using the Calderón formula and an Inverse Fast Fourier Transform.
*/

void ifft_calderon(MatDouble &xi,Array3Complex &qhat)
{
  std::cout << "Computing the conductivity..." << std::endl;

  unsigned ng=qhat.shape()[0];
  unsigned ng2=square(ng);
  unsigned ng3=ng*ng2;
  unsigned i,j,k;
  double nmoon=((double)ng-1.0)/(double)ng;
  double ximax=-xi(0,0);

  //Compute xi^2
  VecDouble normxi2 = vec_norm_2_square_rows_mat(xi);

  //Avoid dividing by zero
  if (ng%2==1)
    {
      normxi2((ng3-1)/2)=1.0;
      qhat[(ng-1)/2][(ng-1)/2][(ng-1)/2]=0.0;
    }

  // Phase
  VecDouble anglephase (ng);
  std::copy(CountItDouble0,CountItDouble0+ng,anglephase.begin());
  anglephase*=-nmoon*Pi;
  VecComplex phase = complexfunvec(VecOp(anglephase,[](double c) -> double { return cos(c); }),VecOp(anglephase,[](double c) -> double { return sin(c); }));

  // Dividing by xi^2, multiplying by 2, dephasing
  for (unsigned p=0;p<ng3;p++)
    {
      i=p/ng2;
      j=(p/ng)%ng; //Same as (p%ng2)/ng
      k=p%ng; //Same as (p%ng2)%ng
      qhat[i][j][k]/=normxi2(p);
      qhat[i][j][k]*=2.0*phase(i)*phase(j)*phase(k);
    }

  // Computing ifft
  fftw_plan plan;
  plan=fftw_plan_dft_3d(ng,ng,ng,reinterpret_cast<fftw_complex *>(&qhat[0][0][0]),reinterpret_cast<fftw_complex *>(&qhat[0][0][0]),FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  // Dephasing, multiplying by coefficient and taking 1-qhat
  dcomplex coeff=pow(nmoon,3)/8.0*dcomplex(cos(3.0*ximax),sin(3.0*ximax));
  for (unsigned p=0;p<ng3;p++)
    {
      i=p/ng2;
      j=(p/ng)%ng; //Same as (p%ng2)/ng
      k=p%ng; //Same as (p%ng2)%ng
      qhat[i][j][k]*=coeff*phase(i)*phase(j)*phase(k);
      qhat[i][j][k]=1.0-qhat[i][j][k];
    }
}
