#ifndef SCHROEDINGER_RHS
#define SCHROEDINGER_RHS

void schroedinger_rhs(VecDouble &X,VecDouble &q,unsigned nnodes,unsigned &ntetrahedra,MatUnsigned &tetrahedra,VecDouble &tetrahedravol,VecUnsigned &dof,VecMatDouble &grad);

#endif // SCHROEDINGER_RHS
