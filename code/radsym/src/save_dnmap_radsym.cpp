#include "radsym.h"

/*
  Saves the Dirichlet-to-Neumann map for a radially symmetric conductivity in dnmapfilename.
*/

void save_dnmap_radsym(unsigned &nd,MatDouble &dnmap,string &dnmapfilename)
{
  std::cout << "Saving the Dirichlet-to-Neumann map in results/dnmap/" << dnmapfilename << std::endl;

  string dnmapfullfilename="results/dnmap/"+dnmapfilename;

  std::ofstream dnmap_file(dnmapfullfilename.c_str(),std::ios::out|std::ios::trunc);

  if (!dnmap_file)
    {
      std::cerr << "Can't open file " << dnmapfullfilename << std::endl;
      exit(1);
    }

  dnmap_file.setf(std::ios::scientific,std::ios::floatfield);
  dnmap_file.precision(15);

  dnmap_file << nd << std::endl;
  for (unsigned i=0;i<dnmap.size1();i++)
    {
      for (unsigned j=0;j<dnmap.size2()-1;j++)
	{
	  dnmap_file << dnmap(i,j) << " ";
	}
      dnmap_file << dnmap(i,dnmap.size2()-1) << std::endl;
    }

  // Closing file
  dnmap_file.close();
}
