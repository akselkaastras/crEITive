#include "eit.h"

/*
  Saves the scattering transform.
*/

void save_q(VecDouble &q,MatDouble &nodes,string &qfilename)
{
  std::cout << "Saving q in results/q/" << qfilename << std::endl;

  string qfullfilename="results/q/"+qfilename;

  std::ofstream q_file(qfullfilename.c_str(),std::ios::out|std::ios::trunc);

  if (!q_file)
    {
      std::cerr << "Can't open file " << qfullfilename << std::endl;
      exit(1);
    }

  q_file.setf(std::ios::scientific,std::ios::floatfield);
  q_file.precision(15);

  unsigned n=nodes.size1();

  q_file << n << std::endl;

  for (unsigned p=0;p<n;p++)
    {
      q_file << nodes(p,0) << " " << nodes(p,1) << " " << nodes(p,2) << " " << q(p) << std::endl;
    }
  // Closing file
  q_file.close();
}
