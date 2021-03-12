#include "eit.h"

/*
  Computes the solution X of the Schrödinger equation.
*/

void solve_schroedinger(CompMatComplexColMaj &M,VecComplex &X,unsigned &nnodes,unsigned &ndof,unsigned &nsurfnodes,VecUnsigned &doftonode,VecUnsigned &surfnodes)
{
  std::cout << "Solving Schrödinger linear system... " << std::endl;

  //Vectors of the compressed matrix
  VecInt Ap (M.index1_data().size());
  copy(M.index1_data().begin(),M.index1_data().end(),Ap.begin());

  VecInt Ai (M.index2_data().size());
  copy(M.index2_data().begin(),M.index2_data().end(),Ai.begin());
  Ai.resize(M.nnz());

  VecComplex Axz (M.value_data().size());
  copy(M.value_data().begin(),M.value_data().end(),Axz.begin());
  Axz.resize(M.nnz());

  VecDouble Ax = real(Axz);
  VecDouble Az = imag(Axz);

  VecDouble Xx = real(X);
  VecDouble Xz = imag(X);

  //umfpack inversion
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  (void) umfpack_zi_symbolic (ndof,ndof,&Ap(0),&Ai(0),&Ax(0),&Az(0),&Symbolic,null,null);
  (void) umfpack_zi_numeric (&Ap(0),&Ai(0),&Ax(0),&Az(0),Symbolic,&Numeric,null,null);
  umfpack_zi_free_symbolic (&Symbolic);
  (void) umfpack_zi_solve (UMFPACK_A,&Ap(0),&Ai(0),&Ax(0),&Az(0),&Xx(0),&Xz(0),&Xx(0),&Xz(0),Numeric,null,null);
  umfpack_zi_free_numeric (&Numeric);

  for (unsigned i=0;i<X.size();i++)
    {
      X(i)=dcomplex(Xx(i),Xz(i));
    }

  VecComplex Xdof=X;
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
