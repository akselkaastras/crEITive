#include "eit.h"

/*
  Computes zeta for a given xi and a given norm sqrt(2)*kappa.

  zeta can be written as
  zeta = zeta_r + i zeta_i
  Where
  zeta_r = -xi/2 + xi^T
  xi^T . xi = 0
  zeta_i . xi = zeta_i . xi^T = 0
  |zeta_r| = |zeta_i| = |zeta| / sqrt(2) = kappa

  zeta is chosen in this way
  Writing xi in spherical coordinates
  xi = |xi| ( sin(theta)cos(phi), sin(theta)sin(phi), cos(theta) )

  zeta_r = -xi/2 + a*( cos(theta)cos(phi), cos(theta)sin(phi), -sin(theta) )
  zeta_i = |zeta_i|*( -sin(phi), cos(phi), 0 )
  With a = sqrt(|zeta_i|^2 - |xi|^2/4) = sqrt(|zeta|^2/2 - |xi|^2/4) = sqrt(kappa^2 - |xi|^2/4)

  Writing zeta = kappa*(k^T + i k)
  kappa >= 0
  k^T,k \in R^3 with k^T.k = 0 and |k|=|k^T|=1

  k^T = zeta_r/kappa = -xi/(2*kappa) + sqrt(kappa^2 - |xi|^2/4)/kappa*( cos(theta)cos(phi), cos(theta)sin(phi), -sin(theta) )
  k = zeta_i/kappa = ( -sin(phi), cos(phi), 0 )

  /!\/!\ /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
  The chosen norm of zeta should verify |zeta|^2 >= |xi|^2/2
  /!\/!\ /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
*/

void zeta_xi(double &pkappa,MatDouble &xi,VecDouble &kappa,MatDouble &k,MatDouble &kT,string &zetamethodname)
{
  double ximax=-xi(0,0);
  //Only xi of norm less than or equal to ximax are considered
  //Find where |xi| <= ximax
  VecDouble vnorm2xi = vec_norm_2_rows_mat(xi);
  VecDouble xinormdiff = ScalVecDouble(kappa.size(),ximax) - vnorm2xi;
  greater_equal_iterator itge (xinormdiff.begin(),greater_equal_0);
  //Done

  if(zetamethodname=="fixed")
    {
      kappa  = pkappa*ScalVecDouble(kappa.size(),ximax/2.0);
    }
  else if (zetamethodname=="proportional")
    {
      kappa=pkappa*vnorm2xi/2.0;
    } 
  //Put to zero when |xi| > ximax
  std::transform(kappa.begin(),kappa.end(),itge,kappa.begin(),std::multiplies<double>());

  k = ZeroMatDouble (k.size1(),k.size2());
  kT = -xi/2.0;

  //Spherical coordinates of xi
  MatDouble xis = cart2sph_rows_mat(xi);
  VecDouble xis0 = boostublas::column(xis,0);
  MatDouble xise = boostublas::subrange(xis,0,xis.size1(),1,xis.size2());
  MatDouble cxise = MatOp(xise,[](double c) -> double { return cos(c); });
  MatDouble sxise = MatOp(xise,[](double c) -> double { return sin(c); });

  //a
  VecDouble a = VecOp(kappa,square<double>)-VecOp(xis0,square<double>)/4.0;
  //Put to zero when |xi| > ximax
  std::transform(a.begin(),a.end(),itge,a.begin(),std::multiplies<double>());
  //Square root
  AutoVecOp(a,[](double c) -> double { return sqrt(c); });


  //k
  boostublas::column(k,0)=-boostublas::column(sxise,1);
  boostublas::column(k,1)=boostublas::column(cxise,1);

  //kT
  //0
  boostublas::column(kT,0)+=boostublas::element_prod(a,boostublas::element_prod(boostublas::column(cxise,0),boostublas::column(cxise,1)));
  //1
  boostublas::column(kT,1)+=boostublas::element_prod(a,boostublas::element_prod(boostublas::column(cxise,0),boostublas::column(sxise,1)));
  //2
  boostublas::column(kT,2)-=boostublas::element_prod(a,boostublas::column(sxise,0));

  //kT of norm 1
  kT = cart2sph_rows_mat(kT);
  //Using former xise cxise sxise for kT now
  xise = boostublas::subrange(kT,0,kT.size1(),1,kT.size2());
  cxise = MatOp(xise,[](double c) -> double { return cos(c); });
  sxise = MatOp(xise,[](double c) -> double { return sin(c); }); 
  boostublas::column(kT,0)=boostublas::element_prod(boostublas::column(sxise,0),boostublas::column(cxise,1));
  boostublas::column(kT,1)=boostublas::element_prod(boostublas::column(sxise,0),boostublas::column(sxise,1));
  boostublas::column(kT,2)=boostublas::column(cxise,0);
}
