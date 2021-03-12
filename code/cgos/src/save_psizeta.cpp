#include "cgos.h"

/*
  Saves the CGO solutions traces.
*/

void save_psizeta(MatComplexColMaj &matpsizeta,double &kappa,VecDouble &akT,VecDouble &theta,VecDouble &phi,VecDouble &tx,VecDouble &phix,string &cgostracesfilename)
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

  unsigned nk=akT.size();
  unsigned nt=theta.size();
  unsigned np=phi.size();
  unsigned ntx=tx.size();
  unsigned nd=ntx-1;
  unsigned nx=2*square(nt);
  unsigned ntm1=nt-1;
  unsigned ngw0=nk*(nt-1)*np;
  unsigned akTk,thetak,phik,kpnk;

  //Header
  cgostraces_file << nd << " " << kappa << " " << nk << " " << nt << " " << np << std::endl;

  //Save quadrature points
  for (unsigned n=0;n<nx-1;n++)
    {
      cgostraces_file << acos(tx(n%ntx)) << " " << phix(n/ntx) << " ";
    }
  cgostraces_file << acos(tx(ntx-1)) << " " << phix(2*ntx-1) << std::endl;

  //Save zeta, psizeta
  //theta=0
  for (unsigned k=0;k<nk;k++)
    {
      cgostraces_file << akT(k) << " " << 0.0 << " " << 0.0 << " ";
      for (unsigned j=0;j<nx-1;j++)
	{
	  cgostraces_file << real(matpsizeta(j,k)) << " " << imag(matpsizeta(j,k)) << " "; 
	}
      cgostraces_file << real(matpsizeta(nx-1,k)) << " " << imag(matpsizeta(nx-1,k)) << std::endl;
    }
  //theta!=0
  for (unsigned k=0;k<ngw0;k++)
    {
      akTk=k%nk; //Same as (k%(nk*ntm1))%nk
      thetak=1+(k/nk)%ntm1; //Same as 1+(k%(nk*ntm1))/nk
      phik=k/(nk*ntm1);
      kpnk=k+nk;
      cgostraces_file << akT(akTk) << " " << theta(thetak) << " " << phi(phik) << " ";
      for (unsigned j=0;j<nx-1;j++)
	{
	  cgostraces_file << real(matpsizeta(j,kpnk)) << " " << imag(matpsizeta(j,kpnk)) << " "; 
	}
      cgostraces_file << real(matpsizeta(nx-1,kpnk)) << " " << imag(matpsizeta(nx-1,kpnk)) << std::endl;
    }

  // Closing file
  cgostraces_file.close();
}
