#ifndef __ETRANS_RND_H
#define __ETRANS_RND_H

#include "config.h"
#include "types.h"

#ifdef USE_MKL
# include <mkl.h>
typedef VSLStreamStatePtr rnd_stream_p;
#else
# include <gsl/gsl_rng.h>
typedef gsl_rng *rnd_stream_p;
#endif

void rnd_uniform_int(int n, int *x, int mx);
void rnd_uniform(int n, real_t *x, real_t mx);
void rnd_gaussian(int n, real_t *x, real_t sigma);
rnd_stream_p rnd_alloc(unsigned int seed);
void rnd_free(rnd_stream_p rnd);

int rnd_init();
void rnd_finish();

#endif /* __ETRANS_RND_H */
