#include "cgos.h"

/*
  Saves the Dirichlet-to-Neumann map for a conductivity equal to 1 in dnmap0filename.
*/

void save_dnmap0(unsigned &nd,MatDoubleColMaj &dnmap0,string &dnmap0filename)
{
  std::cout << "Saving the Dirichlet-to-Neumann map in results/dnmap0/" << dnmap0filename << std::endl;

  string dnmap0fullfilename="results/dnmap0/"+dnmap0filename;

  std::ofstream dnmap0_file(dnmap0fullfilename.c_str(),std::ios::out|std::ios::trunc);

  if (!dnmap0_file)
    {
      std::cerr << "Can't open file " << dnmap0fullfilename << std::endl;
      exit(1);
    }

  dnmap0_file.setf(std::ios::scientific,std::ios::floatfield);
  dnmap0_file.precision(15);

  dnmap0_file << nd << std::endl;
  for (unsigned i=0;i<dnmap0.size1();i++)
    {
      for (unsigned j=0;j<dnmap0.size2()-1;j++)
	{
	  dnmap0_file << dnmap0(i,j) << " ";
	}
      dnmap0_file << dnmap0(i,dnmap0.size2()-1) << std::endl;
    }

  // Closing file
  dnmap0_file.close();
}
