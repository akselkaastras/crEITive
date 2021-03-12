#include "cgos.h"
#include "zeta_rhs.h"
#include "psi_zeta.h"

/*
  Computes the CGO solutions traces for zeta in |zeta| S^2.
*/

void cgos_traces_0(MatDoubleColMaj &diffdtn,MatDoubleColMaj &S0,MatDouble &x,VecDouble &wtrep,double &kappa,MatDouble &k,MatDouble &kT,MatComplexColMaj &matpsizeta)
{
  std::cout << "Computing the CGO solutions traces..." << std::endl;

  int n=diffdtn.size1();
  unsigned ng=k.size1();

  MatDoubleColMaj F0 (n,n);
  MatDoubleColMaj psizeta (n,2);

  // Blas product
  char trans[] = "N";
  double alpha=1.0;
  double beta=0.0;
  dgemm_(trans,trans,n,n,n,alpha,&S0(0,0),n,&diffdtn(0,0),n,beta,&F0(0,0),n);
  //Add identity matrix
  add_identity_matrix(F0);
  // Lapack LU factorization
  VecInt ipiv (n);
  int info;
  dgetrf_(n,n,&F0(0,0),n,&ipiv(0),info);
  ipiv-=ScalVecInt(n,1);
  // For LU inversion
  char side[] = "L";
  char transn[] = "N";
  int m=2;
  char lo[] = "L";
  char diagu[] = "U";
  char up[] = "U";
  char diagn[] = "N";

  std::cout << ng << " points..." << std::endl;

  for (unsigned p=0;p<ng;p++)
    {
      //Compute right hand side of the boundary integral equation (exp(ix.zeta)), stored in psizeta
      zeta_rhs(psizeta,kappa,boostublas::row(k,p),boostublas::row(kT,p),x);
      // Solve the boundary integral equation
      swap_rows(ipiv,psizeta);
      dtrsm_(side,lo,transn,diagu,n,m,alpha,&F0(0,0),n,&psizeta(0,0),n);
      dtrsm_(side,up,transn,diagn,n,m,alpha,&F0(0,0),n,&psizeta(0,0),n);
      //Put in matpsizeta
      boostublas::column(matpsizeta,p)=complexfunmat2(psizeta);
    }
}
