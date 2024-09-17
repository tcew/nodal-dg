#
# Sample system.mk file for gnu compilers
#
CP = cp
MV = mv
CD = cd
MAKE = make

# c++ compiler
CXX = g++
# fortran compiler
FC = gfortran
# loader
LD = g++

# optimization flags passed to all compilers
OPTFLAGS = -O3 -DNDEBUG -Wno-free-nonheap-object
# OPTFLAGS = -g -DDEBUG -D_DEBUG
#OPTFLAGS = -Wall 

# c++ compiler options
CXXOPTIONS = -DUNDERSCORE -fpermissive

# fortran compiler options
FCOPTIONS =

# loader options
LDOPTIONS =

# command to archive the libraries
AR = ar rv

# command to generate an index to the contents of an
# archive
RANLIB = ar -ts
# other possibilities for ranlib
#RANLIB = ranlib
#RANLIB = echo

# indicate that we want to compile the local blas
# and lapack library by setting LOCAL_BLASLAPACK to
# 'libBlasLapack'
# LOCAL_BLASLAPACK = libBlasLapack
# link the local blas and lapack library
BLASLAPACKLIBS = -lBlasLapack -lgfortran
# BLASLAPACKLIBS = -framework Accelerate 
# BLASLAPACKLIBS = -llapack -lblas -lf77blas -latlas
# BLASLAPACKLIBS = -L/Library/Frameworks/Intel_MKL.framework/Libraries/32 -lmkl_lapack32 -lmkl_ia32 

# indicate that we do not want to compile the local blas
# and lapack library by setting LOCAL_BLASLAPACK to ''
#LOCAL_BLASLAPACK =
# link the blas and lapack from ATLAS
#BLASLAPACKLIBS = -L/opt/local/ATLAS/lib/Linux_P4SSE2 \
#                 -llapack -lf77blas -lcblas -latlas

# these are the fortran libs needed to link fortran with
# c++
FORTRANLIBS = -lm
