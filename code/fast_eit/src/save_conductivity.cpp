#include "eit.h"

/*
  Saves the conductivity.
*/

void save_conductivity(MatDouble &nodes,VecDouble &sig,string &conductivityfilename)
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

  conductivity_file << sig.size() << std::endl;
  for (unsigned i=0;i<sig.size();i++)
    {
      conductivity_file << (double) nodes(i,0) << " " << (double) nodes(i,1) << " " << (double) nodes(i,2) << " " << sig(i) << std::endl;
    }

  // Closing file
  conductivity_file.close();
}
