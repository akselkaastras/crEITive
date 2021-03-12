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
#include <boost/numeric/ublas/matrix_sparse.hpp>
//#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/multi_array.hpp>
#include <boost/functional.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/bind.hpp>
#include <fftw3.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <umfpack.h>

namespace boostublas = boost::numeric::ublas;
namespace boostmath = boost::math;

typedef std::string string;

typedef std::complex<double> dcomplex;

typedef boostublas::zero_vector<double> ZeroVecDouble;
typedef boostublas::zero_vector<double> ZeroVecComplex;
typedef boostublas::scalar_vector<int> ScalVecInt;
typedef boostublas::scalar_vector<double> ScalVecDouble;

typedef boostublas::vector<unsigned> VecUnsigned;
typedef boostublas::vector<int> VecInt;
typedef boostublas::vector<double> VecDouble;
typedef boostublas::vector<dcomplex> VecComplex;

typedef boostublas::zero_matrix<double> ZeroMatDouble;
typedef boostublas::zero_matrix<double,boostublas::column_major> ZeroMatDoubleColMaj;
typedef boostublas::zero_matrix<dcomplex,boostublas::column_major> ZeroMatComplexColMaj;

typedef boostublas::matrix<unsigned> MatUnsigned;
typedef boostublas::matrix<double> MatDouble;
typedef boostublas::matrix<double,boostublas::column_major> MatDoubleColMaj;
typedef boostublas::matrix<dcomplex,boostublas::column_major> MatComplexColMaj;

typedef boostublas::vector<boostublas::matrix<double> > VecMatDouble;

typedef boost::multi_array<dcomplex,3> Array3Complex;

//typedef boostublas::permutation_matrix<std::size_t> PermMat;

typedef boostublas::compressed_matrix<dcomplex,boostublas::column_major> CompMatComplexColMaj;

typedef boost::counting_iterator<unsigned> CountItUnsigned;
typedef boost::counting_iterator<double> CountItDouble;

typedef boost::binder2nd<std::greater_equal<double> > binder2nd_greater_equal;

typedef boost::transform_iterator<binder2nd_greater_equal,VecDouble::iterator,boost::use_default,boost::use_default> greater_equal_iterator;

/* #ifndef CONST_IC */
/* #define CONST_IC */
/* const dcomplex Ic (0.0,1.0); */
/* #endif */

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

#ifndef CONST_GREATER_EQUAL_0
#define CONST_GREATER_EQUAL_0
const binder2nd_greater_equal greater_equal_0 = boost::bind2nd(std::greater_equal<double>(),0.0);
#endif

#ifndef EXTERN_DSBEV
#define EXTERN_DSBEV
extern "C" {
  void dsbev_(char *JOBZ, char *UPLO, int &N, int &KD, double *AB, int &LDAB, double *W, double *Z, int &LDZ, double *WORK, int &INFO);
}
#endif

#ifndef EXTERN_ZGESV
#define EXTERN_ZGESV
extern "C" {
  void zgesv_(int &N, int &NRHS, dcomplex *A, int &LDA, int *IPIV, dcomplex *B, int &LDB, int &INFO);
}
#endif

#ifndef EXTERN_ZGEMV
#define EXTERN_ZGEMV
extern "C" {
  void zgemv_(char *TRANS, int &N, int &M, dcomplex &ALPHA, dcomplex *A, int &LDA, dcomplex *X, int &INCX, dcomplex &BETA, dcomplex *Y, int &INCY);
}
#endif

#ifndef EXTERN_DGEMM
#define EXTERN_DGEMM
extern "C" {
  void dgemm_(char *TRANSA, char *TRANSB, int &M, int &N, int &K, double &ALPHA, double *A, int &LDA, double *B, int &LDB, double &BETA, double *C, int &LDC);
}
#endif

#ifndef EXTERN_ZGETRF
#define EXTERN_ZGETRF
extern "C" {
  void zgetrf_(int &M, int &N, dcomplex *A, int &LDA, int *IPIV, int &INFO);
}
#endif

#ifndef EXTERN_ZTRSM
#define EXTERN_ZTRSM
extern "C" {
  void ztrsm_(char *SIDE, char *UPLO, char *TRANSA, char *DIAG, int &M, int &N, dcomplex &ALPHA, dcomplex *A, int &LDA, dcomplex *B, int &LDB);
}
#endif

inline double onetwo(const unsigned &n)
{
  return (n==0) ? 1.0 : 2.0;
}

template<class T>
inline T plus_one(const T &z)
{
  return z+(T)1;
};

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

inline unsigned intsqrt3(const unsigned &n)
{
  if (n<=63)
    {
      if (n==0) return 0;
      else if (n<=7) return 1;
      else if (n<=26) return 2;
      else return 3;
    }
  unsigned sqrt3n=(2*intsqrt(n))/3;
  unsigned sqrt3n2=square(sqrt3n);
  unsigned decsqrt3n;
  while (sqrt3n2>n/sqrt3n)
    {
      decsqrt3n=sqrt3n/3-n/(3*sqrt3n2);
      if (decsqrt3n==0) decsqrt3n=1;
      sqrt3n=sqrt3n-decsqrt3n;
      sqrt3n2=square(sqrt3n);
    }
  return sqrt3n;
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
inline boostublas::matrix<dcomplex,F> complexfunmat(const boostublas::matrix<double,F> &X,const boostublas::matrix<double,F> &Y)
{
  boostublas::matrix<dcomplex,F> Z (X.size1(),X.size2());
  std::transform(X.data().begin(),X.data().end(),Y.data().begin(),Z.data().begin(),complexfun);
  return Z;
};

template<class F>
inline double norm_2_row(const boostublas::matrix<double,F> &M,unsigned i)
{
  return boostublas::norm_2(boostublas::row(M,i));
};

template<class F>
inline VecDouble vec_norm_2_rows_mat(const boostublas::matrix<double,F> &M)
{
  VecDouble V (M.size1());
  std::transform(CountItUnsigned0,CountItUnsigned0+M.size1(),V.begin(),boost::bind1st(norm_2_row<F>,M));
  return V;
};

template<class F>
inline double prec_inner_prod_row(const boostublas::matrix<double,F> &M,unsigned i)
{
  return boostublas::prec_inner_prod(boostublas::row(M,i),boostublas::row(M,i));
};

template<class F>
inline VecDouble vec_norm_2_square_rows_mat(const boostublas::matrix<double,F> &M)
{
  VecDouble V (M.size1());
  std::transform(CountItUnsigned0,CountItUnsigned0+M.size1(),V.begin(),boost::bind1st(prec_inner_prod_row<F>,M));
  return V;
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

template<class F>
inline void prec_inner_prod_mat_mat_row_replace_column(const boostublas::matrix<double,F> &x,const boostublas::matrix<double,F> &M,boostublas::matrix<double,F> &Mx,unsigned j)
{
  boostublas::column(Mx,j)=prec_inner_prod_vec_rows_mat(boostublas::row(x,j),M);
};

template<class F>
inline boostublas::matrix<double,F> prec_inner_prod_mat_mat_rows(const boostublas::matrix<double,F> &x,const boostublas::matrix<double,F> &M)
{
  boostublas::matrix<double,F> Mx (M.size1(),x.size1());
  std::copy(CountItUnsigned0,CountItUnsigned0+x.size1(),boost::make_function_output_iterator(boost::bind(prec_inner_prod_mat_mat_row_replace_column<F>,boost::ref(x),boost::ref(M),boost::ref(Mx),_1)));
  return Mx;
};

inline VecDouble cart2sph(const VecDouble &x)
{
  VecDouble xs = ZeroVecDouble (3);
  xs(0)=sqrt(square(x(0))+square(x(1))+square(x(2)));
  if(xs(0)!=0.0)
    {
      xs(1)=acos(x(2)/xs(0));
      xs(2)=((x(0)==0.0&&x(1)==0.0) ? 0.0 : fmod(atan2(x(1),x(0))+2*Pi,2.0*Pi))
	;
    }
  return xs;
}

template<class F>
inline VecDouble cart2sph_row(const boostublas::matrix<double,F> &M,unsigned i)
{
  return cart2sph(boostublas::row(M,i));
};

template<class F>
inline void cart2sph_replace_row(boostublas::matrix<double,F> &M,unsigned i)
{
  boostublas::row(M,i)=cart2sph(boostublas::row(M,i));
};

template<class F>
inline boostublas::matrix<double,F> cart2sph_rows_mat(const boostublas::matrix<double,F> &M)
{
  boostublas::matrix<double,F> Ms = M;
  std::copy(CountItUnsigned0,CountItUnsigned0+Ms.size1(),boost::make_function_output_iterator(boost::bind1st(cart2sph_replace_row<F>,Ms)));
  return Ms;
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

template<class T,class UnOp>
inline void AutoVecOp(boostublas::vector<T> &X, const UnOp &op)
{
  std::transform(X.begin(),X.end(),X.begin(),op);
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
inline boostublas::vector<T> kron_prod_vec(const boostublas::vector<T> &x,const boostublas::vector<T> &y)
{
  boostublas::vector<T> z (x.size()*y.size());
  boostublas::matrix<T> M = boostublas::outer_prod(x,y);
  std::copy(M.data().begin(),M.data().end(),z.begin());
  return z;
};

template<class It>
inline bool value_iterator(const It &it,unsigned k)
{
  return *(it+k);
};
