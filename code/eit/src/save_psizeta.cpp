#include "eit.h"

/*
  Saves the CGO solutions traces.
*/

void save_psizeta(MatComplexColMaj &matpsizeta,MatDouble &xi,VecDouble &kappa,MatDouble &k,MatDouble &kT,VecDouble &t,VecDouble &phi,string &cgostracesfilename)
{
  std::cout << "Saving the traces of the CGO solutions in results/cgostraces/" << cgostracesfilename << std::endl;

  string cgostracesfullfilename="results/cgostraces/"+cgostracesfilename;

  std::ofstream cgostraces_file(cgostracesfullfilename.c_str(),std::ios::out|std::ios::trunc);

  if (!cgostraces_file)
    {
      std::cerr << "Can't open file " << cgostracesfullfilename << std::endl;
      exit(1);
    }

  cgostraces_file.setf(std::ios::scientific,std::ios::floatfield);
  cgostraces_file.precision(15);

  unsigned nt=t.size();
  unsigned nd=nt-1;
  unsigned nx=2*square(nt);

  double ximax=-xi(0,0);
  unsigned ng3=kappa.size();
  unsigned ng=intsqrt3(ng3);

  //Save quadrature points
  cgostraces_file << nd << " " << ng << " " << ximax << std::endl;
  for (unsigned n=0;n<nx-1;n++)
    {
      cgostraces_file << acos(t(n%nt)) << " " << phi(n/nt) << " ";
    }
  cgostraces_file << acos(t(nt-1)) << " " << phi(2*nt-1) << std::endl;

  //Save xi, zeta, psizeta
  for (unsigned i=0;i<ng3;i++)
    {
      cgostraces_file << xi(i,0) << " " << xi(i,1) << " " << xi(i,2) << " " << kappa(i)*kT(i,0) << " " << kappa(i)*k(i,0) << " " << kappa(i)*kT(i,1) << " " << kappa(i)*k(i,1) << " " << kappa(i)*kT(i,2) << " " << kappa(i)*k(i,2) << " ";
      for (unsigned j=0;j<nx-1;j++)
	{
	  cgostraces_file << real(matpsizeta(j,i)) << " " << imag(matpsizeta(j,i)) << " "; 
	}
      cgostraces_file << real(matpsizeta(nx-1,i)) << " " << imag(matpsizeta(nx-1,i)) << std::endl;
    }

  // Closing file
  cgostraces_file.close();
}
