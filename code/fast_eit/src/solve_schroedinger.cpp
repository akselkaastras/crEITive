#include "eit.h"

/*
  Computes the solution X of the Schrödinger equation.
*/

void solve_schroedinger(CompMatDoubleColMaj &M,VecDouble &X,unsigned &nnodes,unsigned &ndof,unsigned &nsurfnodes,VecUnsigned &doftonode,VecUnsigned &surfnodes)
{
  std::cout << "Solving Schrödinger linear system... " << std::endl;

  //Vectors of the compressed matrix
  VecInt Ap (M.index1_data().size());
  std::copy(M.index1_data().begin(),M.index1_data().end(),Ap.begin());

  VecInt Ai (M.index2_data().size());
  std::copy(M.index2_data().begin(),M.index2_data().end(),Ai.begin());
  Ai.resize(M.nnz());

  VecDouble Ax (M.value_data().size());
  std::copy(M.value_data().begin(),M.value_data().end(),Ax.begin());
  Ax.resize(M.nnz());

  //umfpack inversion
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  (void) umfpack_di_symbolic (ndof,ndof,&Ap(0),&Ai(0),&Ax(0),&Symbolic,null,null);
  (void) umfpack_di_numeric (&Ap(0),&Ai(0),&Ax(0),Symbolic,&Numeric,null,null);
  umfpack_di_free_symbolic (&Symbolic);
  (void) umfpack_di_solve (UMFPACK_A,&Ap(0),&Ai(0),&Ax(0),&X(0),&X(0),Numeric,null,null);
  umfpack_di_free_numeric (&Numeric);

  VecDouble Xdof=X;
  X.resize(nnodes);
  for (unsigned i=0;i<ndof;i++)
    {
      X(doftonode(i))=Xdof(i);
    }
  for (unsigned i=0;i<nsurfnodes;i++)
    {
      X(surfnodes(i))=1.0;
    }
}
