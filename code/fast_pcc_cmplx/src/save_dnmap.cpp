#include "pcc.h"

/*
  Saves the Dirichlet-to-Neumann map in dnmapfilename.
*/

void save_dnmap(unsigned &nd,MatComplex &DNmap,string &dnmapfilename)
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

  dnmap_file << nd << " complex" << std::endl;
  for (unsigned i=0;i<DNmap.size1();i++)
    {
      for (unsigned j=0;j<DNmap.size2()-1;j++)
	{
	  dnmap_file << real(DNmap(i,j)) << " " << imag(DNmap(i,j)) << " ";
	}
      dnmap_file << real(DNmap(i,DNmap.size2()-1)) << " " << imag(DNmap(i,DNmap.size2()-1)) << std::endl;
    }
  dnmap_file.close();
}
