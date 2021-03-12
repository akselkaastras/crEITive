#include "cgos.h"

/*
  Saves the Single Layer Operator Matrix in s0filename.
*/

void save_s0(unsigned &nd,MatDoubleColMaj &S0,string &s0filename)
{
  std::cout << "Saving the Single Layer Operator Matrix in results/s0/" << s0filename << std::endl;

  string s0fullfilename="results/s0/"+s0filename;

  std::ofstream s0_file(s0fullfilename.c_str(),std::ios::out|std::ios::trunc);

  if (!s0_file)
    {
      std::cerr << "Can't open file " << s0fullfilename << std::endl;
      exit(1);
    }

  s0_file.setf(std::ios::scientific,std::ios::floatfield);
  s0_file.precision(15);

  s0_file << nd << std::endl;
  for (unsigned i=0;i<S0.size1();i++)
    {
      for (unsigned j=0;j<S0.size2()-1;j++)
	{
	  s0_file << S0(i,j) << " ";
	}
      s0_file << S0(i,S0.size2()-1) << std::endl;
    }

  // Closing file
  s0_file.close();
}
