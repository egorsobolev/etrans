#include "rnd.h"

#ifdef USE_MKL

#ifdef USE_SINGLE
# define vRngGaussian vsRngGaussian
#else
# define vRngGaussian vdRngGaussian
#endif

void rnd_uniform(rnd_stream_p rnd, int n, int *x, int mx)
{
  viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnd, n, x, 0, mx);
}

void rnd_gaussian(rnd_stream_p rnd, int n, real_t *x, real_t sigma)
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

void rnd_uniform(rnd_stream_p rnd, int n, int *x, int mx)
{
  int i;
  for (i = 0; i < n; ++i)
    x[i] = gsl_rng_uniform_int(rnd, mx);;
}

void rnd_gaussian(rnd_stream_p rnd, int n, real_t *x, real_t sigma)
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
