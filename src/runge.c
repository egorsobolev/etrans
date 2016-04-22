#include "types.h"
#include "chain.h"

#include <rng.h>
#include <math.h>

/* 2o2s1g algorithm */
/* HS Greenside and E Helfand, The Bell System Technical Journal 60 (1981) 1927-1940 */
void rk_2o2s1g(chain_eq_t *c, real_t *wrk, real_t h, real_t *x)
{
  real_t *x0, *g1;
  real_t sqrth;
  int k, n, nn;

  x0 = wrk;
  g1 = x0 + nn;

  sqrth = _sqrt(h) * c->sigF;
  n = c->n;
  nn = 2 * n;

  for (k = 0; k < nn; k++) {
    x0[k] = x[k];
  }
  rng_gaussian(n, x + n, 1.0);
  for (k = n; k < nn; ++k) {
    x[k] = sqrth * x[k] + x0[k];
    g1[k-n] = x[k]; /* du = v */
  }
  c->dv(c, x, g1+n);
  for (k = 0; k < nn; ++k) {
    g1[k] *= h;
    x[k] += 0.5 * g1[k];
    g1[k] += x0[k];
  }
  c->dv(c, g1, x0+n);
  for (k = n; k < nn; ++k) {
    x[k-n] += 0.5 * h * g1[k];
    x[k] += 0.5 * h * x0[k];
  }
}
