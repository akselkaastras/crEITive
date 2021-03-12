#include "radsym.h"

/*
  Saves the radially symmetric conductivity on a radius in the file conductivityfilename.
*/

void save_conductivity(VecDouble &r,VecDouble &conductivity,string &conductivityfilename)
{
  std::cout << "Saving the conductivity on a radius in results/conductivity/" << conductivityfilename << std::endl;

  string conductivityfullfilename="results/conductivity/"+conductivityfilename;

  std::ofstream conductivity_file(conductivityfullfilename.c_str(),std::ios::out|std::ios::trunc);

  if (!conductivity_file)
    {
      std::cerr << "Can't open file " << conductivityfullfilename << std::endl;
      exit(1);
    }

  conductivity_file.setf(std::ios::scientific,std::ios::floatfield);
  conductivity_file.precision(15);

  unsigned n=r.size();

  conductivity_file << "radial" << std::endl;
  conductivity_file << n << std::endl;
  for (unsigned p=0;p<n;p++)
    {
      conductivity_file << r(p) << " " << conductivity(p) << std::endl;
    }

  // Closing file
  conductivity_file.close();
}
