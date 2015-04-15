#include "chain.h"

real_t en_kinetic(int n, const real_t *v)
{
  int k;
  real_t e;
  e = 0.0;
  for (k = 0; k < n; k++)
    e += v[k] * v[k];
  return 0.5 * e;
}

