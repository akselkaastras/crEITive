#include "eit.h"
#include "h_zeta_matrix.h"
#include "zeta_rhs.h"
#include "exp_minus_i_x_xi.h"
#include "psi_zeta.h"
#include "scattering_transform.h"

/*
  Computes the scattering transform on a uniform grid with ngrid points on [-ximax,ximax]^3.
*/

void qhat_zeta(MatDoubleColMaj &diffdtn,MatDoubleColMaj &S0,MatDouble &x,VecDouble &wtrep,MatDouble &xi,const VecDouble &kappa,const MatDouble &k,const MatDouble &kT,MatComplexColMaj &matpsizeta,Array3Complex &qhat)
{
  std::cout << "Computing qhat..." << std::endl;

  int n=diffdtn.size1();
  unsigned ng=qhat.shape()[0];
  unsigned ng2=square(ng);
  unsigned ng3=ng*ng2;

  double ximax=-xi(0,0);

  matpsizeta = ZeroMatComplexColMaj (n,ng3);

  unsigned ii,jj,kk;

  MatDoubleColMaj Hzeta (n,n);
  MatDoubleColMaj Fzeta (n,n);
  MatDoubleColMaj psizeta (n,2);
  VecComplex expmixxipzeta (n);
  dcomplex txizeta;

  // For BLAS product
  char trans[] = "N";
  double alpha=1.0;
  double beta=0.0;

  std::cout << ng3 << " grid points..." << std::endl;

  for (unsigned p=0;p<ng3;p++)
    {
      ii=p/ng2;
      jj=(p/ng)%ng; //Same as (p%ng2)/ng
      kk=p%ng; //Same as (p%ng2)%ng
      qhat[ii][jj][kk]=0.0; //Filled with zeros if not computed
      if (boostublas::norm_2(boostublas::row(xi,p))<=ximax)
	{
	  //Compute Fzeta
	  h_zeta_matrix(Hzeta,kappa(p),boostublas::row(k,p),x,wtrep);
	  // Fzeta=boostublas::prec_prod(S0+Hzeta,diffdtn);
	  // Blas product (faster)
	  Hzeta+=S0;
	  //Blas product
	  dgemm_(trans,trans,n,n,n,alpha,&Hzeta(0,0),n,&diffdtn(0,0),n,beta,&Fzeta(0,0),n);
	  //Add identity matrix
	  add_identity_matrix(Fzeta);
	  //Compute right hand side of the boundary integral equation (exp(ix.zeta)), stored in psizeta
	  zeta_rhs(psizeta,kappa(p),boostublas::row(k,p),boostublas::row(kT,p),x);
	  //exp(-ix.xi)
	  exp_minus_i_x_xi(expmixxipzeta,boostublas::row(xi,p),x);
	  //exp(-ix.(xi+zeta))
	  expmixxipzeta=boostublas::element_div(expmixxipzeta,complexfunmat2(psizeta));
	  //Solve the boundary integral equation
	  psi_zeta(Fzeta,psizeta);
	  //Put in matpsizeta
	  boostublas::column(matpsizeta,p)=complexfunmat2(psizeta);
	  //Scattering transform
	  scattering_transform(txizeta,expmixxipzeta,diffdtn,psizeta,wtrep);
	  //Put in qhat
	  qhat[ii][jj][kk]=txizeta;
	}
    }
}
