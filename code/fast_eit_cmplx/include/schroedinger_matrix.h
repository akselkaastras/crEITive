#ifndef SCHROEDINGER_MATRIX
#define SCHROEDINGER_MATRIX

void schroedinger_matrix(CompMatComplexColMaj &M,VecComplex &q,unsigned nnodes,unsigned &ntetrahedra,MatUnsigned &tetrahedra,VecDouble &tetrahedravol,VecUnsigned &dof,VecMatDouble &grad);

#endif // SCHROEDINGER_MATRIX
