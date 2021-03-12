#include "eit.h"

/*
  Computes the right hand side of the finite elements method for the Schrödinger equation.
  sparse vector but stored as full vector.
*/

void schroedinger_rhs(VecComplex &X,VecComplex &q,unsigned nnodes,unsigned &ntetrahedra,MatUnsigned &tetrahedra,VecDouble &tetrahedravol,VecUnsigned &dof,VecMatDouble &grad)
{
  std::cout << "Computing the Schrödinger right hand side..." << std::endl;

  unsigned dofnti;
  unsigned dofntj;
  VecDouble graddofnti (3);
  VecDouble graddofntj (3);
  X = ZeroVecComplex (X.size());

  for (unsigned nt=0;nt<ntetrahedra;nt++)
    {
      for (unsigned i=0;i<4;i++)
	{
	  dofnti=dof(tetrahedra(nt,i));
	  if (dofnti!=nnodes)
	    {
	      graddofnti=boostublas::row(grad(nt),i);
	      for (unsigned j=0;j<4;j++)
		{
		  dofntj=dof(tetrahedra(nt,j));
		  if (dofntj==nnodes)
		    {
		      graddofntj=boostublas::row(grad(nt),j);
		      X(dofnti)-=(q(nt)/20.0+(double)boostublas::prec_inner_prod(graddofnti,graddofntj))*tetrahedravol(nt);
		    }
		}
	    }
	}
    }
}
