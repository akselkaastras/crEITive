#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <omp.h>
#include <cmath>
#include <complex>
#include <boost/limits.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/multi_array.hpp>

namespace boostublas = boost::numeric::ublas;
namespace boostmath = boost::math;

typedef std::string string;

typedef std::complex<double> dcomplex;

typedef boostublas::vector<double> VecDouble;
typedef boostublas::vector<dcomplex> VecComplex;

typedef boostublas::zero_vector<double> ZeroVecDouble;
typedef boostublas::zero_vector<dcomplex> ZeroVecComplex;

typedef boostublas::scalar_vector<double> ScalVecDouble;

typedef boostublas::matrix<double> MatDouble;
typedef boostublas::matrix<dcomplex> MatComplex;

typedef boostublas::matrix<double,boostublas::column_major> MatDoubleColMaj;
typedef boostublas::matrix<dcomplex,boostublas::column_major> MatComplexColMaj;

typedef boostublas::identity_matrix<double> IdMatDouble;
typedef boostublas::zero_matrix<dcomplex> ZeroMatComplex;

typedef boostublas::zero_matrix<double,boostublas::column_major> ZeroMatDoubleColMaj;
typedef boostublas::zero_matrix<dcomplex,boostublas::column_major> ZeroMatComplexColMaj;

typedef boostublas::vector<MatDouble> VecMatDouble;
typedef boostublas::vector<VecDouble> VecVecDouble;

typedef boost::multi_array<dcomplex,3> Array3Complex;

typedef boostublas::matrix_range<MatComplex> MatRangeComplex;
typedef boostublas::vector_range<VecComplex> VecRangeComplex;

typedef boostublas::permutation_matrix<std::size_t> PermMat;

#ifndef CONDUCTIVITY
#define CONDUCTIVITY

class Conductivity
{
 public:
  Conductivity();
  void SetCenterCoord(unsigned ncoord,double valcoord);
  double GetCenterCoord(unsigned ncoord) const;
  void SetCenter(VecDouble veccenter);
  VecDouble GetCenter() const;
  void SetAxisCoord(unsigned naxis,unsigned ncoord,double valcoord);
  double GetAxisCoord(unsigned naxis,unsigned valcoord) const;
  void SetAxis(unsigned naxis,VecDouble vecaxis);
  VecDouble GetAxis(unsigned naxis) const;
  void SetAxes(MatDouble mataxes);
  MatDouble GetAxes() const;
  void SetRadius(unsigned naxis,double valhalfaxis);
  double GetRadius(unsigned naxis) const;
  void SetRadii(VecDouble vechalfaxes);
  VecDouble GetRadii() const;
  void SetAmplitude(dcomplex valamplitude);
  dcomplex GetAmplitude() const;
  dcomplex GetTwoMu() const;
  void Display() const;
  //
  void SetMaxDegree(unsigned valn);
  unsigned GetMaxDegree() const;
  VecDouble GetElevation() const;
  VecDouble GetCosElevation() const;
  VecDouble GetSinElevation() const;
  VecDouble GetAzimuth() const;
  VecDouble GetCosAzimuth() const;
  VecDouble GetSinAzimuth() const;
  VecDouble GetAlpha() const;
  VecDouble GetBeta() const;
  MatDouble GetSphericalSpherePoints() const;
  MatDouble GetCartesianSpherePoints() const;
  MatDouble GetSurfacePoints() const;
  VecDouble GetJacobian() const;
  MatDouble GetOutwardNormal() const;
  MatDouble GetAssociatedLegendreFunctions() const;
  MatComplex GetComplexExponentials() const;
  MatDouble GetCosElevationRotationsSpherePoints() const;
  MatDouble GetAzimuthRotationsSpherePoints() const;
  MatDouble GetKernelF() const;
  void InitializeData();
 private:
  VecDouble center;
  MatDouble axes;
  VecDouble radii;
  dcomplex amplitude;
  dcomplex twomu;
  //
  unsigned n;
  VecDouble t,ct,st,p,cp,sp;
  VecDouble alpha,beta;
  MatDouble xc,xs,x;
  VecDouble jac;
  MatDouble nu;
  MatDouble lf;
  MatComplex ep;
  MatDouble cttxmoy,ptxmoy,kernf;
};

#endif

typedef boostublas::vector<Conductivity> VecConductivity;

#ifndef CONST_IC
#define CONST_IC
const dcomplex Ic (0.0,1.0);
#endif

#ifndef CONST_EPS
#define CONST_EPS
const double eps = std::numeric_limits<double>::epsilon();
#endif

#ifndef CONST_BIG
#define CONST_BIG
const double big = std::numeric_limits<double>::max();
#endif

#ifndef CONST_PI
#define CONST_PI
const double Pi = boost::math::constants::pi<double>();
#endif

#ifndef CONST_ROOT_PI
#define CONST_ROOT_PI
const double rootPi = boost::math::constants::root_pi<double>();
#endif

#ifndef CONST_ROOT_TWO_PI
#define CONST_ROOT_TWO_PI
const double root2Pi = boost::math::constants::root_two_pi<double>();
#endif

#ifndef CONST_ROOT_FOUR_PI
#define CONST_ROOT_FOUR_PI
const double root4Pi = 2.0*rootPi;
#endif

#ifndef CONST_INV_ROOT_FOUR_PI
#define CONST_INV_ROOT_FOUR_PI
const double invroot4Pi = 0.5/rootPi;
#endif

#ifndef CONST_ROOT_TWO
#define CONST_ROOT_TWO
const double root2 = boost::math::constants::root_two<double>();
#endif

#ifndef EXTERN_DSBEV
#define EXTERN_DSBEV
extern "C" {
  void dsbev_(char *JOBZ, char *UPLO, int &N, int &KD, double *AB, int &LDAB, double *W, double *Z, int &LDZ, double *WORK, int &INFO);
}
#endif

inline unsigned intsqrt(unsigned n)
{
  if (n<=15)
    {
      if (n==0) return 0;
      else if (n<=3) return 1;
      else if (n<=8) return 2;
      else return 3;
    }
  unsigned sqrtn=n/2;
  unsigned decsqrtn;
  while (sqrtn>n/sqrtn)
    {
      decsqrtn=sqrtn/2-n/(2*sqrtn);
      if (decsqrtn==0) decsqrtn=1;
      sqrtn=sqrtn-decsqrtn;
    }
  return sqrtn;
}

inline VecDouble cart2sph(VecDouble x)
{
  VecDouble xs = ZeroVecDouble (3);
  xs(0)=sqrt(pow(x(0),2)+pow(x(1),2)+pow(x(2),2));
  if(xs(0)!=0.0)
    {
      xs(1)=acos(x(2)/xs(0));
      xs(2)=((x(0)==0.0&&x(1)==0.0) ? 0.0 : fmod(atan2(x(1),x(0))+2*Pi,2.0*Pi));
    }
  return xs;
}

inline VecDouble cross_prod(VecDouble x,VecDouble y)
{
  VecDouble z = ZeroVecDouble (3);
  z(0)=x(1)*y(2)-x(2)*y(1);
  z(1)=x(2)*y(0)-x(0)*y(2);
  z(2)=x(0)*y(1)-x(1)*y(0);
  return z;
}

inline double det3(MatDouble M)
{
  return M(0,0)*M(1,1)*M(2,2)+M(0,1)*M(1,2)*M(2,0)+M(0,2)*M(1,0)*M(2,1)-M(0,0)*M(1,2)*M(2,1)-M(0,1)*M(1,0)*M(2,2)-M(0,2)*M(1,1)*M(2,0);
}
