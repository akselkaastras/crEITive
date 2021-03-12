#include "eit.h"
#include "zeta_rhs.h"
#include "exp_minus_i_x_xi.h"
#include "scattering_transform.h"

/*
  Computes the scattering transform with texp approximation on a uniform grid with ngrid points on [-ximax,ximax]^3.
*/

void qhat_zeta_exp(MatDoubleColMaj &diffdtn,MatDouble &x,VecDouble &wtrep,MatDouble &xi,const VecDouble &kappa,const MatDouble &k,const MatDouble &kT,MatComplexColMaj &matpsizeta,Array3Complex &qhat)
{
  std::cout << "Computing qhat..." << std::endl;

  int n=diffdtn.size1();
  unsigned ng=qhat.shape()[0];
  unsigned ng2=square(ng);
  unsigned ng3=ng*ng2;

  double ximax=-xi(0,0);

  matpsizeta = ZeroMatComplexColMaj (n,ng3);

  unsigned ii,jj,kk;

  MatDoubleColMaj expixzeta (n,2);
  VecComplex expmixxipzeta (n);
  dcomplex txizeta;

  std::cout << ng3 << " grid points..." << std::endl;

  for (unsigned p=0;p<ng3;p++)
    {
      ii=p/ng2;
      jj=(p/ng)%ng; //Same as (p%ng2)/ng
      kk=p%ng; //Same as (p%ng2)%ng
      qhat[ii][jj][kk]=0.0; //Filled with zeros if not computed
      if (boostublas::norm_2(boostublas::row(xi,p))<=ximax)
	{
	  //Compute zeta
	  zeta_rhs(expixzeta,kappa(p),boostublas::row(k,p),boostublas::row(kT,p),x);
	  //exp(-ix.xi)
	  exp_minus_i_x_xi(expmixxipzeta,boostublas::row(xi,p),x);
	  //exp(-ix.(xi+zeta))
	  expmixxipzeta=boostublas::element_div(expmixxipzeta,complexfunmat2(expixzeta));
	  //Put exp(ix.zeta) in matpsizeta
	  boostublas::column(matpsizeta,p)=complexfunmat2(expixzeta);
	  //Compute scattering transform
	  scattering_transform(txizeta,expmixxipzeta,diffdtn,expixzeta,wtrep);
	  //Put in 3D array
	  qhat[ii][jj][kk]=txizeta;
	}
    }
}