#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>
#include <cmath>
#include <complex>
#include <boost/limits.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace boostublas = boost::numeric::ublas;

typedef std::string string;

typedef std::complex<double> complex;

typedef boostublas::vector<double> VecDouble;
typedef boostublas::vector<complex> VecComplex;

typedef boostublas::zero_vector<double> ZeroVecDouble;

typedef boostublas::matrix<double> MatDouble;

typedef boostublas::matrix<double,boostublas::column_major> MatDoubleColMaj;

typedef boostublas::zero_matrix<double> ZeroMatDouble;

typedef boostublas::zero_matrix<double,boostublas::column_major> ZeroMatDoubleColMaj;

#ifndef CONST_IC
#define CONST_IC
const complex Ic (0.0,1.0);
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


