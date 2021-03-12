#include "cgos.h"
#include "h_zeta_matrix.h"
#include "zeta_rhs.h"
#include "psi_zeta.h"

/*
  Computes the CGO solutions traces for zeta in |zeta| S^2.
*/

void cgos_traces(MatDoubleColMaj &diffdtn,MatDoubleColMaj &S0,MatDouble &x,VecDouble &wtrep,double &kappa,MatDouble &k,MatDouble &kT,MatComplexColMaj &matpsizeta)
{
  std::cout << "Computing the CGO solutions traces..." << std::endl;

  int n=diffdtn.size1();
  unsigned ng=k.size1();

  MatDoubleColMaj Hzeta (n,n);
  MatDoubleColMaj Fzeta (n,n);
  MatDoubleColMaj psizeta (n,2);

  // For BLAS product
  char trans[] = "N";
  double alpha=1.0;
  double beta=0.0;

  std::cout << ng << " points..." << std::endl;

  for (unsigned p=0;p<ng;p++)
    {
      //Compute Fzeta
      h_zeta_matrix(Hzeta,kappa,boostublas::row(k,p),x,wtrep);
      Hzeta+=S0;
      //Blas product
      dgemm_(trans,trans,n,n,n,alpha,&Hzeta(0,0),n,&diffdtn(0,0),n,beta,&Fzeta(0,0),n);
      //Add identity matrix
      add_identity_matrix(Fzeta);
      //Compute right hand side of the boundary integral equation (exp(ix.zeta)), stored in psizeta
      zeta_rhs(psizeta,kappa,boostublas::row(k,p),boostublas::row(kT,p),x);
      //Solve the boundary integral equation
      psi_zeta(Fzeta,psizeta);
      //Put in matpsizeta
      boostublas::column(matpsizeta,p)=complexfunmat2(psizeta);
    }
}
