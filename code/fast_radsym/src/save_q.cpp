#include "radsym.h"

/*
  Saves the radially symmetric q on a radius in the file qfilename.
*/

void save_q(VecDouble &r,VecDouble &q,string &qfilename)
{
  std::cout << "Saving q on a radius in results/q/" << qfilename << std::endl;

  string qfullfilename="results/q/"+qfilename;

  std::ofstream q_file(qfullfilename.c_str(),std::ios::out|std::ios::trunc);

  if (!q_file)
    {
      std::cerr << "Can't open file " << qfullfilename << std::endl;
      exit(1);
    }

  q_file.setf(std::ios::scientific,std::ios::floatfield);
  q_file.precision(15);

  unsigned n=r.size();

  q_file << "radial" << std::endl;
  q_file << n << std::endl;
  for (unsigned p=0;p<n;p++)
    {
      q_file << r(p) << " " << q(p)  << std::endl;
    }

  // Closing file
  q_file.close();
}
