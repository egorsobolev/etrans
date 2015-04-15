#include "rnd.h"

rnd_stream_p rnd;

#ifdef USE_MKL

#ifdef USE_SINGLE
# define vRngGaussian vsRngGaussian
# define vRngUniform vsRngUniform
#else
# define vRngGaussian vdRngGaussian
# define vRngUniform vdRngUniform
#endif

void rnd_uniform_int(int n, int *x, int mx)
{
  viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnd, n, x, 0, mx);
}

void rnd_uniform(int n, real_t *x, real_t mx)
{
  vRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnd, n, x, 0, mx);
}

void rnd_gaussian(int n, real_t *x, real_t sigma)
{
  vRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, rnd, n, x, 0.0, sigma);
}

rnd_stream_p rnd_alloc(unsigned int seed)
{
  rnd_stream_p rnd;
  if (vslNewStream(&rnd, VSL_BRNG_MT19937, seed) != VSL_STATUS_OK)
    return NULL;

  return rnd;
}

void rnd_free(rnd_stream_p rnd)
{
  vslDeleteStream(&rnd);
}

#else /* USE_MKL */

#include <gsl/gsl_randist.h>

void rnd_uniform_int(int n, int *x, int mx)
{
  int i;
  for (i = 0; i < n; ++i)
    x[i] = gsl_rng_uniform_int(rnd, mx);
}

void rnd_uniform(int n, real_t *x, real_t mx)
{
  int i;
  for (i = 0; i < n; ++i)
    x[i] = gsl_rng_uniform(rnd) * mx;
}

void rnd_gaussian(int n, real_t *x, real_t sigma)
{
  int i;
  for (i = 0; i < n; ++i)
    x[i] = (real_t) gsl_ran_gaussian_ziggurat(rnd, sigma);
}

rnd_stream_p rnd_alloc(unsigned int seed)
{
  rnd_stream_p rnd;
  rnd = gsl_rng_alloc(gsl_rng_mt19937);
  if (!rnd)
    return NULL;
  gsl_rng_set(rnd, seed);
  return rnd;
}

void rnd_free(rnd_stream_p rnd)
{
  gsl_rng_free(rnd);
}

#endif /* USE_MKL */

#ifdef WIN32
# include <process.h>
#define getpid _getpid
typedef int pid_t;

struct timezone
{
  int tz_minuteswest; /* minutes W of Greenwich */
  int tz_dsttime;     /* type of dst correction */
};
struct timeval
{
  time_t tv_sec;
  time_t tv_usec;
};
int gettimeofday(struct timeval *tv, struct timezone *tz);

#else
#include <unistd.h>
#endif

#ifdef MPI
#include <mpi.h>
#endif

int rnd_init()
{
  unsigned int seed;
  struct timeval tm;
  pid_t pid;
#ifdef MPI
  int rank, np;
#endif
  
  pid = getpid();
  gettimeofday(&tm, NULL);
  seed = (int) tm.tv_usec * (int) pid;
#ifdef MPI
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  seed = seed * np + rank;
#endif
  rnd = rnd_alloc(seed);

  return !rnd;
}

void rnd_finish()
{
  rnd_free(rnd);
}
