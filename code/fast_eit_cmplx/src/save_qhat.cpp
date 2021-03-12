#include "eit.h"

/*
  Saves the scattering transform.
*/

void save_qhat(Array3Complex &qhat,MatDouble &xi,string &qhatfilename)
{
  std::cout << "Saving the scattering transform in results/qhat/" << qhatfilename << std::endl;

  string qhatfullfilename="results/qhat/"+qhatfilename;

  std::ofstream qhat_file(qhatfullfilename.c_str(),std::ios::out|std::ios::trunc);

  if (!qhat_file)
    {
      std::cerr << "Can't open file " << qhatfullfilename << std::endl;
      exit(1);
    }

  qhat_file.setf(std::ios::scientific,std::ios::floatfield);
  qhat_file.precision(15);

  double ximax=-xi(0,0);
  unsigned ng=qhat.shape()[0];
  unsigned ng2=square(ng);
  unsigned ng3=ng*ng2;
  unsigned ii,jj,kk;

  qhat_file << ng << " " << ximax << std::endl;
  for (unsigned p=0;p<ng3;p++)
    {
      ii=p/ng2;
      jj=(p/ng)%ng; //Same as (p%ng2)/ng
      kk=p%ng; //Same as (p%ng2)%ng
      qhat_file << xi(p,0) << " " << xi(p,1) << " " << xi(p,2) << " " << real(qhat[ii][jj][kk]) << " " << imag(qhat[ii][jj][kk]) << std::endl;
    }
  // Closing file
  qhat_file.close();
}
