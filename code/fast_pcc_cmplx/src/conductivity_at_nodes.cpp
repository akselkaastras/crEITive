#include "pcc.h"

/*
  Computes the values of the conductivity at given nodes.
*/

void conductivity_at_nodes(VecComplex &sigmanodes,MatDouble &nodes,VecConductivity &sigma)
{
  std::cout << "Computing the conductivity at the given nodes ..." << std::endl;
  std::cout << "(not used to compute the Dirichlet-to-Neumann map)" << std::endl;
  std::cout << "(intended to be saved for possible comparisons)" << std::endl;

  unsigned nnodes=nodes.size1();
  unsigned nsigma=sigma.size();
  sigmanodes = ZeroVecComplex (nnodes);

  //For conductivity parameters
  VecVecDouble vcenter (nsigma);
  VecMatDouble maxes (nsigma);
  VecVecDouble vradii (nsigma);
  VecComplex valamp (nsigma);

  // For computations
  double islessone;
  VecDouble nodesloop (3);
  VecDouble nodesloop2 (3);
  VecVecDouble vradii2 (nsigma);

  for (unsigned k=0;k<nsigma;k++)
    {
      vcenter(k)=sigma(k).GetCenter();
      maxes(k)=sigma(k).GetAxes();
      vradii(k)=sigma(k).GetRadii();
      valamp(k)=sigma(k).GetAmplitude();

      vradii2(k)=vradii(k);
      vradii2(k)(0)=pow(vradii2(k)(0),2);
      vradii2(k)(1)=pow(vradii2(k)(1),2);
      vradii2(k)(2)=pow(vradii2(k)(2),2);
    }

  for (unsigned inodes=0;inodes<nnodes;inodes++)
    {
      for (unsigned k=0;k<nsigma;k++)
	{
	  nodesloop=boostublas::row(nodes,inodes);
	  nodesloop-=vcenter(k); //Translation
	  nodesloop=boostublas::prec_prod(maxes(k),nodesloop); //Rotation
	  nodesloop2(0)=pow(nodesloop(0),2);
	  nodesloop2(1)=pow(nodesloop(1),2);
	  nodesloop2(2)=pow(nodesloop(2),2);

	  islessone=nodesloop2(0)/vradii2(k)(0)+nodesloop2(1)/vradii2(k)(1)+nodesloop2(2)/vradii2(k)(2);

	  if (islessone<=1.0)
	    {
	      sigmanodes(inodes)+=(valamp(k)-1.0);
	    }
	}
      sigmanodes(inodes)+=1.0;
    }
}
