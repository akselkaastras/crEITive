#include "eit.h"

/*
  Reads a dnmap file.
*/

void read_dnmap(string &filename,unsigned &nd,MatDoubleColMaj &dnmap)
{

  string fullfilename="results/dnmap/"+filename;

  // Opening file
  std::ifstream dnmap_file(fullfilename.c_str(),std::ios::in);

  if (!dnmap_file)
    {
      std::cerr << "Can't open file " << fullfilename << std::endl;
      exit(1);
    }

  std::cout << "File " << fullfilename << " successfully opened" << std::endl;


  std::cout << "Reading the Dirichlet-to-Neumann map file..." << std::endl;

  string dnmap_line;
  //Reading dnmap header (maximal degree of spherical harmonics)
  dnmap_file >> nd;
  getline(dnmap_file,dnmap_line);

  std::cout << "Maximal degree of spherical harmonics:" << std::endl;
  std::cout << "nd=" << nd << std::endl;

  unsigned n=2*square(nd+1);
  dnmap.resize(n,n);

  //Reading dnmap values
  for (unsigned i=0;i<n;i++)
    {
      for (unsigned j=0;j<n;j++)
	{
	  dnmap_file >> dnmap(i,j);
	}
      getline(dnmap_file,dnmap_line);
    }

  // Closing file
  dnmap_file.close();
}
