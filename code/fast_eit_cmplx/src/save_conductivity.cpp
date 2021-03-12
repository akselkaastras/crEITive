#include "eit.h"

/*
  Saves the conductivity.
*/

void save_conductivity(MatDouble &nodes,VecComplex &sig,string &conductivityfilename)
{
  std::cout << "Saving the conductivity in results/conductivity/" << conductivityfilename << std::endl;

  string conductivityfullfilename="results/conductivity/"+conductivityfilename;

  std::ofstream conductivity_file(conductivityfullfilename.c_str(),std::ios::out|std::ios::trunc);

  if (!conductivity_file)
    {
      std::cerr << "Can't open file" << conductivityfullfilename << std::endl;
      exit(1);
    }

  conductivity_file.setf(std::ios::scientific,std::ios::floatfield);
  conductivity_file.precision(15);

  conductivity_file << sig.size() << " complex" << std::endl;
  for (unsigned i=0;i<sig.size();i++)
    {
      conductivity_file << nodes(i,0) << " " << nodes(i,1) << " " << nodes(i,2) << " " << real(sig(i)) << " " << imag(sig(i)) << std::endl;
    }

  // Closing file
  conductivity_file.close();
}
