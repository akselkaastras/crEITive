#include "eit.h"
#include "zeta_rhs.h"
#include "exp_minus_i_x_xi.h"
#include "psi_zeta.h"
#include "scattering_transform.h"

/*
  Computes the scattering transform with t0 approximation on a uniform grid with ngrid points on [-ximax,ximax]^3.
*/

void qhat_zeta_0(MatComplexColMaj &diffdtn,MatDoubleColMaj &rdiffdtn,MatDoubleColMaj &idiffdtn,MatDoubleColMaj &S0,MatDouble &x,VecDouble &wtrep,MatDouble &xi,const VecDouble &kappa,const MatDouble &k,const MatDouble &kT,MatComplexColMaj &matpsizeta,Array3Complex &qhat)
{
  std::cout << "Computing qhat..." << std::endl;

  int n=diffdtn.size1();
  unsigned ng=qhat.shape()[0];
  unsigned ng2=square(ng);
  unsigned ng3=ng*ng2;

  double ximax=-xi(0,0);

  matpsizeta = ZeroMatComplexColMaj (n,ng3);

  unsigned ii,jj,kk;

  MatDoubleColMaj rF0 (n,n);
  MatDoubleColMaj iF0 (n,n);
  VecComplex psizeta (n);
  VecComplex expmixxipzeta (n);
  dcomplex txizeta;

  // F0=boostublas::prec_prod(S0,diffdtn);
  // Blas product (faster)
  char trans[] = "N";
  double alpha=1.0;
  double beta=0.0;
  //Blas product
  dgemm_(trans,trans,n,n,n,alpha,&S0(0,0),n,&rdiffdtn(0,0),n,beta,&rF0(0,0),n);
  dgemm_(trans,trans,n,n,n,alpha,&S0(0,0),n,&idiffdtn(0,0),n,beta,&iF0(0,0),n);
  //Add identity matrix
  add_identity_matrix(rF0);
  // LU factorization
  // PermMat pm (n);
  // int info=lu_factorize(F0,pm);
  // Lapack LU factorization (faster)
  VecInt ipiv (n);
  int info;
  MatComplexColMaj F0 = complexfunmat(rF0,iF0);
  zgetrf_(n,n,&F0(0,0),n,&ipiv(0),info);
  ipiv-=ScalVecInt(n,1);
  // For LU inversion
  dcomplex alphac = 1.0;
  char side[] = "L";
  char transn[] = "N";
  int m=1;
  char lo[] = "L";
  char diagu[] = "U";
  char up[] = "U";
  char diagn[] = "N";

  std::cout << ng3 << " grid points..." << std::endl;

  for (unsigned p=0;p<ng3;p++)
    {
      ii=p/ng2;
      jj=(p/ng)%ng; //Same as (p%ng2)/ng
      kk=p%ng; //Same as (p%ng2)%ng
      qhat[ii][jj][kk]=0.0; //Filled with zeros if not computed
      if (boostublas::norm_2(boostublas::row(xi,p))<=ximax)
	{
	  //Compute right hand side of the boundary integral equation (exp(ix.zeta)), stored in psizeta
	  zeta_rhs(psizeta,kappa(p),boostublas::row(k,p),boostublas::row(kT,p),x);
	  //exp(-ix.xi)
	  exp_minus_i_x_xi(expmixxipzeta,boostublas::row(xi,p),x);
	  //exp(-ix.(xi+zeta))
	  expmixxipzeta=boostublas::element_div(expmixxipzeta,psizeta);
	  // Solve the boundary integral equation
	  // swap_rows(pm,psizeta); //Permute rows of right hand side
	  // psizeta=boostublas::solve(F0,boostublas::solve(F0,psizeta,boostublas::unit_lower_tag()),boostublas::upper_tag()); //Solve system
	  // Solve the boundary integral equation using LAPACK (faster)
	  //Permute rows of right hand side
	  swap_rows(ipiv,psizeta);
	  ztrsm_(side,lo,transn,diagu,n,m,alphac,&F0(0,0),n,&psizeta(0),n);
	  ztrsm_(side,up,transn,diagn,n,m,alphac,&F0(0,0),n,&psizeta(0),n);
	  //Put in matpsizeta
	  boostublas::column(matpsizeta,p)=psizeta;
	  //Compute scattering transform
	  scattering_transform(txizeta,expmixxipzeta,diffdtn,psizeta,wtrep);
	  //Put in 3D array
	  qhat[ii][jj][kk]=txizeta;
	}
    }
}
