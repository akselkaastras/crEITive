#include "pcc.h"
#include "gauss_legendre.h"
#include "phi_values.h"
#include "associated_legendre_functions.h"
#include "exp_i_m_phi.h"

Conductivity::Conductivity()
{
  center = ZeroVecDouble (3);
  axes = IdMatDouble (3);
  radii = ScalVecDouble (3,1.0);
  amplitude = 1.0;
  n=0;
}
void Conductivity::SetCenterCoord(unsigned ncoord,double valcoord)
{
  center(ncoord)=valcoord;
}
double Conductivity::GetCenterCoord(unsigned ncoord) const
{
  return center(ncoord);
}
void Conductivity::SetCenter(VecDouble veccenter)
{
  center=veccenter;
}
VecDouble Conductivity::GetCenter() const
{
  return center;
}
void Conductivity::SetAxisCoord(unsigned naxis,unsigned ncoord,double valcoord)
{
  axes(naxis,ncoord)=valcoord;
}
double Conductivity::GetAxisCoord(unsigned naxis,unsigned ncoord) const
{
  return axes(naxis,ncoord);
}
void Conductivity::SetAxis(unsigned naxis,VecDouble vecaxis)
{
  boostublas::row(axes,naxis)=vecaxis;
}
VecDouble Conductivity::GetAxis(unsigned naxis) const
{
  return boostublas::row(axes,naxis);
}
void Conductivity::SetAxes(MatDouble mataxes)
{
  axes=mataxes;
}
MatDouble Conductivity::GetAxes() const
{
  return axes;
}
void Conductivity::SetRadius(unsigned naxis,double valradius)
{
  radii(naxis)=valradius;
}
double Conductivity::GetRadius(unsigned naxis) const
{
  return radii(naxis);
}
void Conductivity::SetRadii(VecDouble vecradii)
{
  radii=vecradii;
}
VecDouble Conductivity::GetRadii() const
{
  return radii;
}
void Conductivity::SetAmplitude(double valamplitude)
{
  amplitude=valamplitude;
  twomu=2.0*(1.0-amplitude)/(1.0+amplitude);
}
double Conductivity::GetAmplitude() const
{
  return amplitude;
}

double Conductivity::GetTwoMu() const
{
  return twomu;
}
void Conductivity::Display() const
{
  if (radii(0)==radii(1)&&radii(1)==radii(2))
    {
      std::cout << "Ball center:" << std::endl;
      std::cout << center << std::endl;
      std::cout << "Ball radius:" << std::endl;
      std::cout << radii(0) << std::endl;
    }
  else
    {
      std::cout << "Ellipsoid center:" << std::endl;
      std::cout << center << std::endl;
      std::cout << "Ellipsoid axes direction vectors:" << std::endl;
      std::cout << axes << std::endl;
      std::cout << "Ellipsoid radii:" << std::endl;
      std::cout << radii << std::endl;
    }
  std::cout << "Amplitude:" << std::endl;
  std::cout << amplitude << std::endl;
  std::cout << "Number of quadrature points:" << std::endl;
  std::cout << 2*(n+1)*(n+1) << " (n=" << n << ")"<< std::endl;
}

void Conductivity::SetMaxDegree(unsigned valn)
{
  n=valn;
}
unsigned Conductivity::GetMaxDegree() const
{
  return n;
}
VecDouble Conductivity::GetElevation() const
{
  return t;
}
VecDouble Conductivity::GetCosElevation() const
{
  return ct;
}
VecDouble Conductivity::GetSinElevation() const
{
  return st;
}
VecDouble Conductivity::GetAzimuth() const
{
  return p;
}
VecDouble Conductivity::GetCosAzimuth() const
{
  return cp;
}
VecDouble Conductivity::GetSinAzimuth() const
{
  return sp;
}
VecDouble Conductivity::GetAlpha() const
{
  return alpha;
}
VecDouble Conductivity::GetBeta() const
{
  return beta;
}
MatDouble Conductivity::GetSphericalSpherePoints() const
{
  return xs;
}
MatDouble Conductivity::GetCartesianSpherePoints() const
{
  return xc;
}
MatDouble Conductivity::GetSurfacePoints() const
{
  return x;
}
VecDouble Conductivity::GetJacobian() const
{
  return jac;
}
MatDouble Conductivity::GetOutwardNormal() const
{
  return nu;
}
MatDouble Conductivity::GetAssociatedLegendreFunctions() const
{
  return lf;
}
MatComplex Conductivity::GetComplexExponentials() const
{
  return ep;
}
MatDouble Conductivity::GetCosElevationRotationsSpherePoints() const
{
  return cttxmoy;
}
MatDouble Conductivity::GetAzimuthRotationsSpherePoints() const
{
  return ptxmoy;
}
MatDouble Conductivity::GetKernelF() const
{
  return kernf;
}
void Conductivity::InitializeData()
{
  t.resize(n+1);
  ct.resize(n+1);
  st.resize(n+1);
  p.resize(2*(n+1));
  cp.resize(2*(n+1));
  sp.resize(2*(n+1));
  alpha.resize(n+1);
  beta.resize(n+1);
  xs.resize(2*(n+1)*(n+1),2);
  xc.resize(2*(n+1)*(n+1),3);
  //
  x.resize(2*(n+1)*(n+1),3);
  jac.resize(2*(n+1)*(n+1));
  nu.resize(2*(n+1)*(n+1),3);
  //
  lf.resize(n+1,((n+1)*(n+2))/2);
  ep.resize(2*(n+1),2*n+1);
  //
  cttxmoy.resize(2*(n+1)*(n+1),2*(n+1)*(n+1));
  ptxmoy.resize(2*(n+1)*(n+1),2*(n+1)*(n+1));
  kernf.resize(2*(n+1)*(n+1),2*(n+1)*(n+1));
  //
  gauss_legendre(ct,alpha);
  
  alpha*=Pi/((double)(n+1));
  std::transform(ct.begin(),ct.end(),t.begin(),[](double c) -> double { return acos(c); });
  std::transform(t.begin(),t.end(),st.begin(),[](double c) -> double { return sin(c); });
  //
  beta = ZeroVecDouble (n+1);
  double plp1,pl,plm1;
  for (unsigned it=0;it<n+1;it++)
    {
      plm1=boostmath::legendre_p(0,ct(it));
      pl=boostmath::legendre_p(1,ct(it));
      beta(it)=plm1+pl;
      for (unsigned l=1;l<n;l++)
	{
	  plp1=boostmath::legendre_next(l,ct(it),pl,plm1);
	  beta(it)+=plp1;
	  plm1=pl;
	  pl=plp1;
	}
      beta(it)*=alpha(it);
    }
  //
  phi_values(p);
  std::transform(p.begin(),p.end(),cp.begin(),[](double c) -> double { return cos(c); });
  std::transform(p.begin(),p.end(),sp.begin(),[](double c) -> double { return sin(c); });
  //
  if (det3(axes)<0) boostublas::row(axes,2)=-boostublas::row(axes,2);
  //
  unsigned ni;
  for (unsigned ip=0;ip<2*(n+1);ip++)
    {
      for (unsigned it=0;it<n+1;it++)
	{
	  ni=it+(n+1)*ip;
	  xs(ni,0)=t(it);xs(ni,1)=p(ip);
	  xc(ni,0)=st(it)*cp(ip);xc(ni,1)=st(it)*sp(ip);xc(ni,2)=ct(it);
	  x(ni,0)=radii(0)*xc(ni,0);x(ni,1)=radii(1)*xc(ni,1);x(ni,2)=radii(2)*xc(ni,2);
	  boostublas::row(x,ni)=center+boostublas::prec_prod(boostublas::trans(axes),boostublas::row(x,ni));
	  boostublas::row(nu,ni)=xc(ni,0)*radii(1)*radii(2)*boostublas::row(axes,0)+xc(ni,1)*radii(2)*radii(0)*boostublas::row(axes,1)+xc(ni,2)*radii(0)*radii(1)*boostublas::row(axes,2);
	  jac(ni)=boostublas::norm_2(boostublas::row(nu,ni));
	  boostublas::row(nu,ni)/=jac(ni);
	}
    }
  //
  associated_legendre_functions(lf,ct,n);
  //
  exp_i_m_phi(ep,p);
  //
  VecDouble txmoy (3);
  VecDouble txmoys (3);
  VecDouble y (3);
  VecDouble xloop (3);
  double jactxmoy,normdiffxquad;
  VecDouble diffx (3);
  VecDouble cp2 (2*(n+1));
  VecDouble sp2 (2*(n+1));
  for (unsigned jp=0;jp<2*(n+1);jp++)
    {
      cp2(jp)=pow(cp(jp),2);
      sp2(jp)=pow(sp(jp),2);
    }
  unsigned nj;
  for (unsigned ip=0;ip<2*(n+1);ip++) 
    {
      for (unsigned it=0;it<n+1;it++) 
	{
	  ni=it+(n+1)*ip;
	  y=boostublas::row(xc,ni);
	  for (unsigned jp=0;jp<2*(n+1);jp++)
	    {
	      for (unsigned jt=0;jt<n+1;jt++) 
		{
		  nj=jt+(n+1)*jp;
		  y=boostublas::row(xc,ni);
		  txmoy(0)=(cp2(jp)*ct(jt)+sp2(jp))*y(0)+cp(jp)*sp(jp)*(ct(jt)-1.0)*y(1)+cp(jp)*st(jt)*y(2);
		  txmoy(1)=cp(jp)*sp(jp)*(ct(jt)-1.0)*y(0)+(sp2(jp)*ct(jt)+cp2(jp))*y(1)+sp(jp)*st(jt)*y(2);
		  txmoy(2)=-cp(jp)*st(jt)*y(0)-sp(jp)*st(jt)*y(1)+ct(jt)*y(2);
		  txmoys=cart2sph(txmoy);
		  cttxmoy(ni,nj)=cos(txmoys(1));
		  ptxmoy(ni,nj)=txmoys(2);

		  jactxmoy=boostublas::norm_2(txmoy(0)*radii(1)*radii(2)*boostublas::row(axes,0)+txmoy(1)*radii(2)*radii(0)*boostublas::row(axes,1)+txmoy(2)*radii(0)*radii(1)*boostublas::row(axes,2));

		  xloop=boostublas::row(xc,nj)-txmoy;
		  normdiffxquad=norm_2(xloop);
		  xloop(0)*=radii(0);
		  xloop(1)*=radii(1);
		  xloop(2)*=radii(2);
		  diffx=boostublas::prec_prod(boostublas::trans(axes),xloop);
		  diffx/=pow(boostublas::norm_2(diffx),3);

		  kernf(ni,nj)=-alpha(jt)*beta(it)*jactxmoy*normdiffxquad*boostublas::prec_inner_prod(diffx,boostublas::row(nu,nj));

		}
	    }
	}
    }
}
