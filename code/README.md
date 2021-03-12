Compilation of code
===================

The code requires a set of libraries:

- BLAS/LAPACK
  These should be provided through non-NetLib
  libraries for best performance.

  Typically OpenBLAS is a good candidate.

  Usually just doing:

    module load openblas

  would be fine.

- GSL is a suite of routines for calling
  scientific specific methods.

  GSL is located here /appl/gsl/2.5/

- fftw3 it is located here in /lib and need not be loaded.
  It is not compatible with the required compiler.

- suitesparse (containing UMFPACK)

   module load suitesparse


Generally the required libraries should be
edited in arch.make. This is *imported* for all
sub-codes and thus one can tune details in one place.


So to make it compile, simply do:

$> module load openblas suitesparse
$> make

   
