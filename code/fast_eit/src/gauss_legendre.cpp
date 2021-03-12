#include "eit.h"

/*
  Computes the nodes and weights for the Gauss-Legendre quadrature formula
  of order n on [-1,1] with the Golub-Welsch algorithm. That is, writing
  the recursion formula on the Legendre monic polynomials Q_n :
  For m>=1, Q_{m+1} + (B_m - X)Q_m + A_m Q_{m-1} = 0

  The nodes are given by the eigenvalues of the Jacobi tridigonal matrix :

  (   B_0     sqrt(A_1)     0          ...             0       )
  ( sqrt(A_1)   B_1      sqrt(A_2)     ...             0       )
  (    0      sqrt(A_2)    B_3         ...             0       )
  (   ...       ...        ...         ...       sqrt(A_{n-1}) )
  (    0         0         ...     sqrt(A_{n-1})    B_{n-1})   )

  Considering v_p, 1<=p<=n, normalized eigenvectors of the matrix, the
  weights are given by :
  mu_0 v_p(1)^2
  Where mu_0 is the integral of the weight function \omega for the
  orthogonality of the polynomials.

  In our particular case of Legendre polynomials :
  (m+1) P_{m+1} - (2m+1) X P_m + m P_{m-1} = 0
  We have :
  A_m = m^2/(4m^2-1)
  B_m = 0
  \omega = 1
  mu_0 = 2

  x is the vector of nodes.
  w is the vector of weights.
*/

void gauss_legendre(VecDouble &x,VecDouble &w)
{
  std::cout << "Computing the Gauss-Legendre nodes and weights..." << std::endl;

  /*
    Matrix m with two lines to store the diagonals of
    the real symmetric band Jacobi matrix.
    The first row starting to second column is the upper diagonal.
    The second row is the main diagonal (zeros).
  */

  MatDoubleColMaj m = ZeroMatDoubleColMaj (2,x.size());
  for (unsigned j=1;j<x.size();j++)
    {
      m(0,j)=(double)j/sqrt(4.0*(double)square(j)-1.0);
    }

  //Lapack
  char jobz[] = "V"; //Computes eigenvalues and eigenvectors
  char uplo[] = "U"; //Upper triangle of the Jacobi matrix is stored
  int n = x.size(); //Size of the Jacobi matrix
  int kd = 1; //Number of superdiagonals
  int ldab = 2; //Leading dimension of m
  MatDoubleColMaj z (x.size(),x.size()); //Orthonormal eigenvectors
  int ldz = x.size(); //Leading dimension of z
  VecDouble work (3*x.size()-2); //Workspace
  int info; //Information integer
  dsbev_(jobz,uplo,n,kd,&m(0,0),ldab,&x(0),&z(0,0),ldz,&work(0),info);

  //Weights
  for (unsigned j=0;j<w.size();j++)
    {
      w(j)=2.0*pow(z(0,j),2.0);
    }

  //Symmetrize
  for (unsigned j=0;j<x.size();j++)
    {
      x(j)=(x(j)-x(x.size()-j-1))/2.0;
    }
  for (unsigned j=0;j<w.size();j++)
    {
      w(j)=(w(j)+w(w.size()-j-1))/2.0;
    }
}
