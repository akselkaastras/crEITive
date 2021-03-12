#include "radsym.h"

/*
  Saves the eigenvalues of the Dirichlet-to-Neumann map for a radially symmetric conductivity in dnmapeigenvaluesfilename, as well as the estimated error.
*/

void save_dnmap_eigenvalues(double &stepsize,VecDouble &eigenvalues,VecDouble &eigenvalueserror,string &dnmapeigenvaluesfilename)
{
  std::cout << "Saving the eigenvalues of the Dirichlet-to-Neumann map in results/dnmapeigenvalues/" << dnmapeigenvaluesfilename << std::endl;

  string dnmapeigenvaluesfullfilename="results/dnmapeigenvalues/"+dnmapeigenvaluesfilename;

  std::ofstream dnmapeigenvalues_file(dnmapeigenvaluesfullfilename.c_str(),std::ios::out|std::ios::trunc);

  if (!dnmapeigenvalues_file)
    {
      std::cerr << "Can't open file " << dnmapeigenvaluesfullfilename << std::endl;
      exit(1);
    }

  dnmapeigenvalues_file.setf(std::ios::scientific,std::ios::floatfield);
  dnmapeigenvalues_file.precision(15);

  unsigned nd=eigenvalues.size()-1;

  dnmapeigenvalues_file << nd << " " << stepsize << std::endl;
  for (unsigned p=0;p<=nd;p++)
    {
      dnmapeigenvalues_file << eigenvalues(p) << " " << eigenvalueserror(p) << std::endl;
    }

  // Closing file
  dnmapeigenvalues_file.close();
}
