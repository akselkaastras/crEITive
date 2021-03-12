#include "pcc.h"

/*
  Saves the conductivity.
*/

void save_conductivity(MatDouble &nodes,VecDouble &sigmanodes,string &conductivityfilename)
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

  conductivity_file << sigmanodes.size() << std::endl;
  for (unsigned i=0;i<sigmanodes.size();i++)
    {
      conductivity_file << nodes(i,0) << " " << nodes(i,1) << " " << nodes(i,2) << " " << sigmanodes(i) << std::endl;
    }

  // Closing file
  conductivity_file.close();
}
