#include "cgos.h"

/*
  Reads a s0 file.
*/

void read_s0(string &filename,MatDoubleColMaj &S0)
{

  string fullfilename="results/s0/"+filename;

  // Opening file
  std::ifstream s0_file(fullfilename.c_str(),std::ios::in);

  if (!s0_file)
    {
      std::cerr << "Can't open file " << fullfilename << std::endl;
      exit(1);
    }

  std::cout << "File " << fullfilename << " successfully opened" << std::endl;


  std::cout << "Reading the Single Layer Operator file..." << std::endl;

  string s0_line;
  unsigned nd;
  //Reading s0 header (maximal degree of spherical harmonics)
  s0_file >> nd;
  getline(s0_file,s0_line);
  unsigned n=2*square(nd+1);
  if (n!=S0.size1())
    {
      std::cerr << "Bad size for the Single Layer Operator" << std::endl;
      exit(1);
    }

  //Reading S0 values
  for (unsigned i=0;i<n;i++)
    {
      for (unsigned j=0;j<n;j++)
	{
	  s0_file >> S0(i,j);
	}
      getline(s0_file,s0_line);
    }

  // Closing file
  s0_file.close();
}
