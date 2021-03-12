#include "eit.h"

/*
  Computes the matrix of the finite elements method for the Schrödinger equation.
  Sparse matrix but stored as full matrix.
*/

void schroedinger_matrix(CompMatDoubleColMaj &M,VecDouble &q,unsigned nnodes,unsigned &ntetrahedra,MatUnsigned &tetrahedra,VecDouble &tetrahedravol,VecUnsigned &dof,VecMatDouble &grad)
{
  std::cout << "Computing the Schrödinger matrix..." << std::endl;

  unsigned dofnti;
  unsigned dofntj;
  VecDouble graddofnti (3);
  VecDouble graddofntj (3);

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
		  if (dofntj!=nnodes)
		    {
		      graddofntj=boostublas::row(grad(nt),j);
		      if (dofnti==dofntj)
			{
			  M(dofnti,dofntj)+=(q(nt)/10.0+boostublas::prec_inner_prod(graddofnti,graddofntj))*tetrahedravol(nt);
			}
		      else
			{
			  M(dofnti,dofntj)+=(q(nt)/20.0+boostublas::prec_inner_prod(graddofnti,graddofntj))*tetrahedravol(nt);
			}
		    }
		}
	    }
	}
    }
  //Symmetrize M to decrease error
  M=(M+boostublas::trans(M))/2.0;
}
