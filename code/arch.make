

# Variables that should be defined are:
#  CXX      -- C++ compiler
#  CXXFLAGS -- C++ flags
#  LIBS     -- Linker flags
#
# This program requires the following libraries which should be linkable through
# the LIBS variable:
#  Boost
#  gsl
#  SuiteSparse (Umfpack and AMD)
#  FFTW
#  BLAS + LAPACK
#  OpenMP


# # Define compiler
# CXX = clang++
# CXXFLAGS = -O3 -Wall $(INC) -fopenmp

# # GSL
# GSL_PATH = /appl/gsl/2.5
# INC += -I$(GSL_PATH)/include
# LIBS += -L$(GSL_PATH)/lib -Wl,-rpath=$(GSL_PATH)/lib -lgsl -lgslcblas

# # Suitesparse
# SUITESPARSE_PATH = /appl/SuiteSparse/5.1.2-sl73
# LIBS += -L$(SUITESPARSE_PATH)/lib -Wl,-rpath=$(SUITESPARSE_PATH)/lib -lumfpack -lamd

# # FFTW
# #FFTW_PATH = /appl/fftw/3.3.8/${CPUTYPEV}
# # If you have a specific path for FFTW, please adapt
# #LIBS += -L$(FFTW_PATH)/lib -Wl,-rpath=$(FFTW_PATH)/lib -lfftw3
# LIBS += -lfftw3

# # BLAS + LAPACK
# LIBS += -lopenblas


# setting ostype 
OSTYPE    = $(shell uname)


ifeq ($(OSTYPE), Linux)
##### linux #####

	#GSL_PATH = /appl/gsl/2.5
	#SUITESPARSE_PATH = /appl/SuiteSparse/5.1.2-sl73


	#INC += -I/usr/local
	#INC += -I$(GSL_PATH)/include

	#LIBS += -L$(GSL_PATH)/lib -Wl,-rpath=$(GSL_PATH)/lib -lgsl -lgslcblas
	#LIBS += -L$(SUITESPARSE_PATH)/lib -Wl,-rpath=$(SUITESPARSE_PATH)/lib -lumfpack -lamd
	#LIBS += -lfftw3
	#LIBS += -L$(OPENBLAS_PATH)/lib -Wl,-rpath=$(OPENBLAS_PATH)/lib -lopenblas

	GSL_PATH = /appl/gsl/2.5
	SUITESPARSE_PATH = /appl/SuiteSparse/5.1.2-sl73
	#FFTW_PATH = /appl/fftw/3.3.8/${CPUTYPEV}

	INC += -I$(GSL_PATH)/include
	INC += -I$(SUITESPARSE_PATH)/include

	LDFLAGS += $(addprefix -L, $(subst :, ,${LD_LIBRARY_PATH}))
	LDFLAGS += -L$(GSL_PATH)/lib -Wl,-rpath=$(GSL_PATH)/lib -lgsl -lgslcblas
	#LDFLAGS += -L$(FFTW_PATH)/lib -Wl,-rpath=$(FFTW_PATH)/lib -lfftw3
	LDFLAGS += -lfftw3
	LDFLAGS += -L$(SUITESPARSE_PATH)/lib -Wl,-rpath=$(SUITESPARSE_PATH)/lib -lumfpack -lamd
	LDFLAGS += -lopenblas



	CXX = g++
	CXXFLAGS = -O3 -Wall $(INC)


else ifeq ($(OSTYPE), Darwin)
##### macOS #####

	GSL_PATH = /usr/local/Cellar/gsl/2.6
	SUITESPARSE_PATH = /usr/local/Cellar/suite-sparse/5.8.1
	OPENBLAS_PATH = /usr/local/Cellar/openblas/0.3.13
	FFTW_PATH = /usr/local/Cellar/fftw/3.3.9

	INC += -I$(GSL_PATH)/include
	INC += -I/usr/local

	LDFLAGS += $(addprefix -L, $(subst :, ,${LD_LIBRARY_PATH}))
	LDFLAGS += -L$(GSL_PATH)/lib -Wl -lgsl -lgslcblas
	LDFLAGS += -lfftw3
	LDFLAGS += -L$(SUITESPARSE_PATH)/lib -Wl -lumfpack -lamd
	LDFLAGS += -L$(OPENBLAS_PATH)/lib  -Wl -lopenblas

	CXX = clang++
	CXXFLAGS = -std=c++11 -DARMA_DONT_PRINT_CXX11_WARNING -O3 -Wall $(INC) -Xclang -fopenmp

endif


