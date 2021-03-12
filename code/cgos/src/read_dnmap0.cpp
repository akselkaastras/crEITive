#include "cgos.h"

/*
  Reads a dnmap file.
*/

void read_dnmap0(string &filename,MatDoubleColMaj &dnmap0)
{

  string fullfilename="results/dnmap0/"+filename;

  // Opening file
  std::ifstream dnmap0_file(fullfilename.c_str(),std::ios::in);

  if (!dnmap0_file)
    {
      std::cerr << "Can't open file " << fullfilename << std::endl;
      exit(1);
    }

  std::cout << "File " << fullfilename << " successfully opened" << std::endl;


  std::cout << "Reading the Dirichlet-to-Neumann map file for a conductivity equal to 1..." << std::endl;

  string dnmap0_line;
  unsigned nd;
  //Reading dnmap0 header (maximal degree of spherical harmonics)
  dnmap0_file >> nd;
  getline(dnmap0_file,dnmap0_line);
  unsigned n=2*square(nd+1);
  if (n!=dnmap0.size1())
    {
      std::cerr << "Bad size for the Dirichlet-to-Neumann map" << std::endl;
      exit(1);
    }

  //Reading dnmap0 values
  for (unsigned i=0;i<n;i++)
    {
      for (unsigned j=0;j<n;j++)
	{
	  dnmap0_file >> dnmap0(i,j);
	}
      getline(dnmap0_file,dnmap0_line);
    }

  // Closing file
  dnmap0_file.close();
}
