#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <complex>
#include <functional>
#include <algorithm>
#include <sys/time.h>
#include <sys/resource.h>
#include <boost/limits.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/functional.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/bind.hpp>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

namespace boostublas = boost::numeric::ublas;
namespace boostmath = boost::math;

typedef std::string string;

typedef std::complex<double> dcomplex;

typedef boostublas::zero_vector<double> ZeroVecDouble;
typedef boostublas::scalar_vector<int> ScalVecInt;
typedef boostublas::scalar_vector<double> ScalVecDouble;

typedef boostublas::vector<int> VecInt;
typedef boostublas::vector<double> VecDouble;
typedef boostublas::vector<dcomplex> VecComplex;

typedef boostublas::zero_matrix<double,boostublas::column_major> ZeroMatDoubleColMaj;

typedef boostublas::matrix<unsigned> MatUnsigned;
typedef boostublas::matrix<double> MatDouble;
typedef boostublas::matrix<double,boostublas::column_major> MatDoubleColMaj;
typedef boostublas::matrix<dcomplex,boostublas::column_major> MatComplexColMaj;

typedef boost::counting_iterator<unsigned> CountItUnsigned;
typedef boost::counting_iterator<double> CountItDouble;

#ifndef CONST_EPS
#define CONST_EPS
const double eps = std::numeric_limits<double>::epsilon();
#endif

#ifndef CONST_PI
#define CONST_PI
const double Pi = boost::math::constants::pi<double>();
#endif

#ifndef CONST_ROOT_PI
#define CONST_ROOT_PI
const double rootPi = boost::math::constants::root_pi<double>();
#endif

#ifndef CONST_INV_ROOT_FOUR_PI
#define CONST_INV_ROOT_FOUR_PI
const double invroot4Pi = 0.5/rootPi;
#endif

#ifndef CONST_COUNTIT_UNSIGNED_0
#define CONST_COUNTIT_UNSIGNED_0
const CountItUnsigned CountItUnsigned0 = CountItUnsigned(0);
#endif

#ifndef CONST_COUNTIT_DOUBLE_0
#define CONST_COUNTIT_DOUBLE_0
const CountItDouble CountItDouble0 = CountItDouble(0);
#endif

#ifndef CONST_COUNTIT_DOUBLE_1
#define CONST_COUNTIT_DOUBLE_1
const CountItDouble CountItDouble1 = CountItDouble(1);
#endif

#ifndef EXTERN_DSBEV
#define EXTERN_DSBEV
extern "C" {
  void dsbev_(char *JOBZ, char *UPLO, int &N, int &KD, double *AB, int &LDAB, double *W, double *Z, int &LDZ, double *WORK, int &INFO);
}
#endif

#ifndef EXTERN_DGESV
#define EXTERN_DGESV
extern "C" {
  void dgesv_(int &N, int &NRHS, double *A, int &LDA, int *IPIV, double *B, int &LDB, int &INFO);
}
#endif

#ifndef EXTERN_DGEMM
#define EXTERN_DGEMM
extern "C" {
  void dgemm_(char *TRANSA, char *TRANSB, int &M, int &N, int &K, double &ALPHA, double *A, int &LDA, double *B, int &LDB, double &BETA, double *C, int &LDC);
}
#endif

#ifndef EXTERN_DGETRF
#define EXTERN_DGETRF
extern "C" {
  void dgetrf_(int &M, int &N, double *A, int &LDA, int *IPIV, int &INFO);
}
#endif

#ifndef EXTERN_DTRSM
#define EXTERN_DTRSM
extern "C" {
  void dtrsm_(char *SIDE, char *UPLO, char *TRANSA, char *DIAG, int &M, int &N, double &ALPHA, double *A, int &LDA, double *B, int &LDB);
}
#endif

inline double onetwo(const unsigned &n)
{
  return (n==0) ? 1.0 : 2.0;
}

template<class T>
inline T square(const T &z)
{
  return z*z;
};

inline unsigned intsqrt(const unsigned &n)
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

inline dcomplex complexfun(const double &x,const double &y)
{
  return dcomplex(x,y);
}

template<class T,class F>
  inline void add_identity_matrix(boostublas::matrix<T,F> &M)
{
  for (unsigned md=0;md<M.size1();md++)
    {
      M(md,md)+=1.0;
    }
};

inline VecComplex complexfunvec(const VecDouble &X,const VecDouble &Y)
{
  VecComplex Z (X.size());
  std::transform(X.begin(),X.end(),Y.begin(),Z.begin(),complexfun);
  return Z;
}

template<class F>
inline VecComplex complexfunmat2(const boostublas::matrix<double,F> &XY)
{
  return complexfunvec(boostublas::column(XY,0),boostublas::column(XY,1));
};

template<class F>
inline double prec_inner_prod_vec_row(const VecDouble &x,const boostublas::matrix<double,F> &M,unsigned i)
{
  return boostublas::prec_inner_prod(x,boostublas::row(M,i));
};

template<class F>
inline VecDouble prec_inner_prod_vec_rows_mat(const VecDouble &x,const boostublas::matrix<double,F> &M)
{
  VecDouble V (M.size1());
  std::transform(CountItUnsigned0,CountItUnsigned0+M.size1(),V.begin(),boost::bind(prec_inner_prod_vec_row<F>,boost::ref(x),boost::ref(M),_1));
  return V;
};

template<class T,class UnOp>
  inline boostublas::vector<T> VecOp(const boostublas::vector<T> &X, const UnOp &op)
{
  boostublas::vector<T> OpX (X.size());
  std::transform(X.begin(),X.end(),OpX.begin(),op);
  return OpX;
};

template<class T,class F,class UnOp>
  inline boostublas::matrix<T,F> MatOp(const boostublas::matrix<T,F> &M, const UnOp &op)
{
  boostublas::matrix<T,F> OpM (M.size1(),M.size2());
  std::transform(M.data().begin(),M.data().end(),OpM.data().begin(),op);
  return OpM;
};

template<class T,class F,class UnOp>
  inline void AutoMatOp(boostublas::matrix<T,F> &M, const UnOp &op)
{
  std::transform(M.data().begin(),M.data().end(),M.data().begin(),op);
};

template<class T>
inline void replace_range(boostublas::vector<T> &xrep,const boostublas::vector<T> &x,unsigned i)
{
  boostublas::subrange(xrep,i*x.size(),(i+1)*x.size())=x;
};

template<class T>
inline boostublas::vector<T> rep_vec_1(const boostublas::vector<T> &x,unsigned n)
{
  boostublas::vector<T> xrep (n*x.size());
  std::copy(CountItUnsigned0,CountItUnsigned0+n,boost::make_function_output_iterator(boost::bind(replace_range<T>,boost::ref(xrep),boost::ref(x),_1)));
  return xrep;
};

template<class T>
inline void replace_slice(boostublas::vector<T> &xrep,const boostublas::vector<T> &x,unsigned n,unsigned i)
{
  boostublas::subslice(xrep,i,n,x.size())=x;
};

template<class T>
inline boostublas::vector<T> rep_vec_2(const boostublas::vector<T> &x,unsigned n)
{
  boostublas::vector<T> xrep (n*x.size());
  std::copy(CountItUnsigned0,CountItUnsigned0+n,boost::make_function_output_iterator(boost::bind(replace_slice<T>,boost::ref(xrep),boost::ref(x),boost::ref(n),_1)));
  return xrep;
};

template<class T>
inline void replace_range_rep_vec_range(boostublas::vector<T> &xrep,const boostublas::vector<T> &x,unsigned nel,unsigned nrep,unsigned i)
{
  replace_range(xrep,rep_vec_2((VecDouble)boostublas::subrange(x,nel*i,nel*(i+1)),nrep),i);
};

template<class T>
inline boostublas::vector<T> rep_vec(const boostublas::vector<T> &x,unsigned nel,unsigned nrep)
{
  boostublas::vector<T> xrep (nrep*x.size());
  std::copy(CountItUnsigned0,CountItUnsigned0+x.size()/nel,boost::make_function_output_iterator(boost::bind(replace_range_rep_vec_range<T>,boost::ref(xrep),boost::ref(x),boost::ref(nel),boost::ref(nrep),_1)));
  return xrep;
};

template<class T>
inline boostublas::vector<T> kron_prod_vec(const boostublas::vector<T> &x,const boostublas::vector<T> &y)
{
  boostublas::vector<T> z (x.size()*y.size());
  boostublas::matrix<T> M = boostublas::outer_prod(x,y);
  std::copy(M.data().begin(),M.data().end(),z.begin());
  return z;
};
