#include "eit.h"
#include "zeta_rhs.h"
#include "exp_minus_i_x_xi.h"
#include "psi_zeta.h"
#include "scattering_transform.h"
#include "butter_coef.h"

/*
  Computes the scattering transform with t0 approximation on a uniform grid with ngrid points on [-ximax,ximax]^3.
*/

void qhat_zeta_0(MatDoubleColMaj &diffdtn,MatDoubleColMaj &S0,MatDouble &x,VecDouble &wtrep,MatDouble &xi,const VecDouble &kappa,const MatDouble &k,const MatDouble &kT,MatComplexColMaj &matpsizeta,Array3Complex &qhat, unsigned &N, unsigned &A)
{
  std::cout << "Computing qhat..." << std::endl;

  int n=diffdtn.size1();
  unsigned ng=qhat.shape()[0];
  unsigned ng2=square(ng);
  unsigned ng3=ng*ng2;

  double ximax=-xi(0,0);

  // Butterworth low-pass filter
  double coef = 1;
  double ximax2 = 0;
  double steplen = 0;
  double steplen2 = 0;
  double ximaxsteplen = 0;
  double omegac = 0;

  if (N > 0 && A > 0) {
    ximax2 = ximax*ximax;
    steplen = Pi*(double)(ng-1)/((double)ng);
    steplen2 = steplen*steplen;
    ximaxsteplen = ximax*steplen;
    omegac = pow((double)(A*A-1),(double) 1/N)/(ximax2); // constant for computing Butterworth coefficient
}  else if (N == 0){
}  else {
    std::cout << "Choose A>0 and N>0 or just N = 0" << std::endl;
    exit(1);
}

  matpsizeta = ZeroMatComplexColMaj (n,ng3);

  MatDoubleColMaj F0 (n,n);

  // F0=boostublas::prec_prod(S0,diffdtn);
  // Blas product (faster)
  char trans[] = "N";
  double alpha=1.0;
  double beta=0.0;
  //Blas product
  dgemm_(trans,trans,n,n,n,alpha,&S0(0,0),n,&diffdtn(0,0),n,beta,&F0(0,0),n);
  //Add identity matrix
  add_identity_matrix(F0);
  // LU factorization
  // PermMat pm (n);
  // int info=lu_factorize(F0,pm);
  // Lapack LU factorization (faster)
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

  std::cout << ng3 << " grid points..." << std::endl;

#pragma omp parallel for default(none) shared(diffdtn,x,wtrep,xi,kappa,k,kT,matpsizeta,qhat,n,ng,ng2,ng3,ximax,F0,alpha,ipiv,side,transn,m,lo,diagu,up,diagn,ximax2,steplen,steplen2,ximaxsteplen,N,omegac,coef)
  for (unsigned p=0;p<ng3;p++)
    {
      unsigned ii=p/ng2;
      unsigned jj=(p/ng)%ng; //Same as (p%ng2)/ng
      unsigned kk=p%ng; //Same as (p%ng2)%ng

      MatDoubleColMaj psizeta (n,2);
      VecComplex expmixxipzeta (n);
      dcomplex txizeta;

      qhat[ii][jj][kk]=0.0; //Filled with zeros if not computed
      if (boostublas::norm_2(boostublas::row(xi,p))<=ximax)
	{
	  //Compute right hand side of the boundary integral equation (exp(ix.zeta)), stored in psizeta
	  zeta_rhs(psizeta,kappa(p),boostublas::row(k,p),boostublas::row(kT,p),x);
	  //exp(-ix.xi)
	  exp_minus_i_x_xi(expmixxipzeta,boostublas::row(xi,p),x);
	  //exp(-ix.(xi+zeta))
	  expmixxipzeta=boostublas::element_div(expmixxipzeta,complexfunmat2(psizeta));
	  // Solve the boundary integral equation
	  // swap_rows(pm,psizeta); //Permute rows of right hand side
	  // psizeta=boostublas::solve(F0,boostublas::solve(F0,psizeta,boostublas::unit_lower_tag()),boostublas::upper_tag()); //Solve system
	  // Solve the boundary integral equation using LAPACK (faster)
	  //Permute rows of right hand side
	  swap_rows(ipiv,psizeta);
	  dtrsm_(side,lo,transn,diagu,n,m,alpha,&F0(0,0),n,&psizeta(0,0),n);
	  dtrsm_(side,up,transn,diagn,n,m,alpha,&F0(0,0),n,&psizeta(0,0),n);
	  //Put in matpsizeta
	  boostublas::column(matpsizeta,p)=complexfunmat2(psizeta);
	  //Compute scattering transform
	  scattering_transform(txizeta,expmixxipzeta,diffdtn,psizeta,wtrep);
    // Butter coefficient multiplication
    if (N != 0) {
      coef = butter_coef(ii,jj,kk,ximax,ximax2,steplen,steplen2,ximaxsteplen,N,omegac);
    }
	  //Put in 3D array
	  qhat[ii][jj][kk]=coef*txizeta;
	}
    }
}