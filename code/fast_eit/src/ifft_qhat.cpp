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
  double ximax=-xi(0,0);
  double nmoon1=(2*(double)ximax)/((double)ng-1.0);
  double nmoon = ((double)ng-1.0)/(double)ng;
  double ximaxo = Pi*((double)ng-1)*nmoon/2.0;
  double scale = (double) ximax/((double) ximaxo);
  //std::cout << ximax << "  is truncrad" << std::endl;

  //std::cout << qhat[3][4][5] << "  is qhat before" << std::endl;
  //std::cout << qhat[2][3][5] << "  is qhat before" << std::endl;
  // Phase
  VecDouble anglephase (ng);
  VecDouble anglephase1 (ng);
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
  //std::cout << qhat[3][4][5] << "  is qhat" << std::endl;
  //std::cout << qhat[2][3][5] << "  is qhat" << std::endl;
  // Computing ifft
  fftw_plan plan;
  plan=fftw_plan_dft_3d(ng,ng,ng,reinterpret_cast<fftw_complex *>(&qhat[0][0][0]),reinterpret_cast<fftw_complex *>(&qhat[0][0][0]),FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
   // std::cout << qhat[3][4][5] << "  is qhat after fft" << std::endl;
  //std::cout << qhat[2][3][5] << "  is qhat after fft" << std::endl;
  // Dephasing and multiplying by coefficient
  //dcomplex coeff=pow((double)nmoon/((double) Pi),3)/8.0*dcomplex(cos(3.0*ximax),sin(3.0*ximax));
  dcomplex coeff = pow(nmoon*scale,3)/8.0*dcomplex(cos(3.0*ximaxo),sin(3.0*ximaxo));
  for (unsigned p=0;p<ng3;p++)
    {
      i=p/ng2;
      j=(p/ng)%ng; //Same as (p%ng2)/ng
      k=p%ng; //Same as (p%ng2)%ng
      qhat[i][j][k]*=coeff*phase(i)*phase(j)*phase(k);
    }

  //    std::cout << qhat[3][4][5] << "  is qhat final" << std::endl;
  //std::cout << qhat[2][3][5] << "  is qhat final" << std::endl;
}
