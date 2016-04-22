#ifndef __ETRANS_TYPES_H
#define __ETRNAS_TYPES_H

#include "config.h"

#ifdef USE_SINGLE
typedef float real_t;
# define MPI_ET_REAL MPI_FLOAT

# define _sqrt(X) sqrtf(X)
# define _sin(X) sinf(X)
# define _cos(X) cosf(X)
# define _fabs(X) fabsf(X)
# define _log(X) logf(X)
# define _exp(X) expf(X)
# define _ceil(X) ceilf(X)
# define _fmod(X) fmodf(X)
# define cblas_copy cblas_scopy
# define cblas_scal cblas_sscal
# define cblas_dot cblas_sdot
# define cblas_nrm2 cblas_snrm2
# define cblas_axpy cblas_saxpy
# define clapack_lamch clapack_slamch
# define clapack_stevx clapack_sstevx
# define thomas1 thomas1s

# define rng_gaussian rng_gaussian_s
# define rng_uniform rng_uniform_s
#else
typedef double real_t;
# define MPI_ET_REAL MPI_DOUBLE

# define _sqrt(X) sqrt(X)
# define _sin(X) sin(X)
# define _cos(X) cos(X)
# define _fabs(X) fabs(X)
# define _log(X) log(X)
# define _exp(X) exp(X)
# define _ceil(X) ceil(X)
# define _fmod(X) fmodf(X)
# define cblas_copy cblas_dcopy
# define cblas_scal cblas_dscal
# define cblas_dot cblas_ddot
# define cblas_nrm2 cblas_dnrm2
# define cblas_axpy cblas_daxpy
# define clapack_lamch clapack_dlamch
# define clapack_stevx clapack_dstevx
# define thomas1 thomas1d

# define rng_gaussian rng_gaussian_d
# define rng_uniform rng_uniform_d
#endif

#endif /* __ETRANS_TYPES_H */
