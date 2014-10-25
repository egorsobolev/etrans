#include "config.h"
#include "osc.h"

#include <time.h>
#include <string.h>
#include <math.h>

#ifdef USE_MKL
# include <mkl.h>
#else
# include <gsl/gsl_cblas.h>
#endif

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

#include "etrans.h"

void f(const chain_eq_t *chain, const real_t *x, real_t *dx)
{
  real_t *dv;
  const real_t *v;
  int i;
  const oscilator_t *o = (oscilator_t *) chain;

  dv = dx + o->n;
  v = x + o->n;
  cblas_copy(o->n, v, 1, dx, 1);
  dv[0] = -o->upr[0] * x[0] - o->tren[0] * v[0] + o->xi[0] * x[1];
  for (i = 1; i < o->n1; ++i)
    dv[i] = -o->upr[i] * x[i] - o->tren[i] * v[i] + o->xi[i] * (x[i + 1] + x[i - 1]);
  dv[o->n1] = -o->upr[o->n1] * x[o->n1] - o->tren[o->n1] * v[o->n1] + o->xi[o->n1] * x[o->n1 - 1];
}

void f2(const chain_eq_t *chain, const real_t *x, real_t *dx, real_t *b2)
{
  real_t *dv;
  const real_t *v;
  int i;
  const oscilator_t *o = (oscilator_t *) chain;

  dv = dx + o->n;
  v = x + o->n;
  cblas_copy(o->n, v, 1, dx, 1);
  dv[0] = -o->upr[0] * x[0] - o->tren[0] * v[0] + o->xi[0] * x[1] - o->lambda[0] * b2[0];
  for (i = 1; i < o->n1; ++i)
    dv[i] = -o->upr[i] * x[i] - o->tren[i] * v[i] + o->xi[i] * (x[i + 1] + x[i - 1]) - o->lambda[i] * b2[i];
  dv[o->n1] = -o->upr[o->n1] * x[o->n1] - o->tren[o->n1] * v[o->n1] + o->xi[o->n1] * x[o->n1 - 1] - o->lambda[o->n1] * b2[o->n1];
}

void runge(chain_eq_t *chain, real_t h, real_t *x)
{
  real_t sqrth, h_2, *y, *k1, *k2;
  int i, m;
  oscilator_t *o = (oscilator_t *) chain;

  sqrth = _sqrt(h);
  h_2 = 0.5 * h;
  m = 2 * o->n;
  y = o->runge_temp;
  k1 = y + m;
  k2 = k1 + m;

  cblas_copy(m, x, 1, y, 1);

  rnd_gaussian(o->rstr, o->n, x + o->n, 1.0);
  /* memset(x + o->n, 0, o->n * sizeof(float)); */

  for (i = o->n; i < m; ++i)
    x[i] = sqrth * o->sig[i - o->n] * x[i] + y[i];
  f(chain, x, k1);
  cblas_axpy(m, h, k1, 1, y, 1);
  f(chain, y, k2);
  cblas_axpy(m, h_2, k1, 1, x, 1);
  cblas_axpy(m, h_2, k2, 1, x, 1);
}

void runge2(chain_eq_t *chain, real_t h, real_t *x, real_t *b0, real_t *b1)
{
  real_t sqrth, h_2, *y, *k1, *k2;
  int i, m;
  oscilator_t *o = (oscilator_t *) chain;

  sqrth = _sqrt(h);
  h_2 = 0.5 * h;
  m = 2 * o->n;
  y = o->runge_temp;
  k1 = y + m;
  k2 = k1 + m;

  cblas_copy(m, x, 1, y, 1);

  rnd_gaussian(o->rstr, o->n, x + o->n, 1.0);
  /*memset(x + o->n, 0, o->n * sizeof(float));*/

  for (i = o->n; i < m; ++i)
    x[i] = sqrth * o->sig[i - o->n] * x[i] + y[i];
  f2(chain, x, k1, b0);
  cblas_axpy(m, h, k1, 1, y, 1);
  f2(chain, y, k2, b1);
  cblas_axpy(m, h_2, k1, 1, x, 1);
  cblas_axpy(m, h_2, k2, 1, x, 1);
}

void osc_equilibrate(chain_eq_t *chain, real_t h, real_t *x, int nstep, int rstep)
{
  int i, n;
  oscilator_t *o = (oscilator_t *) chain;

  if (rstep > 0) {
    rnd_uniform(o->rstr, 1, &n, rstep);
    n += nstep;
  } else
    n = nstep;
  for (i = 0; i < n; ++i)
    runge(chain, h, x);
}
/*
void osc_integrate(oscilator_t *osc, real_t h, real_t *x, int nskip, int nout, real_t *u, real_t *v)
{
  int i, j, k;
  k = 0;
  for (j = 0; j < nout; ++j) {
    for (i = 0; i < nskip; ++i)
      runge(osc, h, x);
    cblas_copy(osc->n, x, 1, u + k, 1);
    cblas_copy(osc->n, x + osc->n, 1, v + k, 1);
    k += osc->n;
  }
}
*/
void osc_x0_rand(chain_eq_t *chain, real_t *x)
{
  int i, n2;
  double a, b, xi, *u;
  oscilator_t *o = (oscilator_t *) chain;

  n2 = 2 * o->n;
  if (o->kt == 0.0 || o->K == 0.0) {
    memset(x, 0, n2 * sizeof(real_t));
    return;
  }
  o->sv = _sqrt(0.5 * o->kt);
  rnd_gaussian(o->rstr, o->n, x + o->n, o->sv);

  /*  printf("K=%g D=%g\n", o->K, o->D);*/

  if (o->D == 0.0) {
    o->su = _sqrt(0.5 * o->kt / o->K);
    b = 1.0;
    a = 0.0;
  } else {
    b = 2.0 * o->D;
    a = _sqrt(o->K * (o->K + 2.0 * b)); 

    o->su = _sqrt(0.5 * o->kt / a);
    xi = (o->K + b - a) / b;

    /*    printf("sigma=%g xi=%g\n", sigma, xi);*/

    a = _sqrt(xi);
    b = _sqrt(1.0 - xi);
  }
  rnd_gaussian(o->rstr, o->n, x, o->su);

  if (o->D == 0.0)
    return;

  for (i = 1; i < o->n; ++i) {
    x[i] *= b;
    x[i] += x[i-1] * a; 
  }
  /*
  for (i = 0; i < o->n; ++i) {
    if (o->K == 0.0)
      x[i] = x[i + o->n] = 0.0;
    else if (o->D == 0.0)
      x[i] /= _sqrt(o->K);
    else
      x[i] *= (1.0 / a)
  }
  */
}
void osc_free(chain_eq_t *chain)
{
  oscilator_t *o = (oscilator_t *) chain;
  rnd_free(o->rstr);
  free(o->xi);
  free(chain);
}


oscilator_t *osc_init(int n, real_t temp, real_t tren, real_t upr, real_t xi, real_t lambda)
{
  static real_t e0 = 0.261838952;
  static real_t t0 = 100.0;
  int i, a;
  unsigned int seed;
  struct timeval tm;
  pid_t pid;
#ifdef MPI
  int rank, np;
#endif
  oscilator_t *o;

  o = (oscilator_t *) malloc(sizeof(oscilator_t));
  if (!o)
    return NULL;

  o->K = upr;
  o->D = xi;
  o->L = tren;
  o->chi = lambda;

  o->kt = e0 / t0 * temp;
  o->n = n;
  o->n1 = n - 1;
  o->xi = (real_t *) malloc(11 * n * sizeof(real_t));
  if (!o->xi) {
    free(o);
    return NULL;
  }

  o->period = 2.0 * M_PI / _sqrt(upr);

  pid = getpid();
  gettimeofday(&tm, NULL);
  seed = (int) tm.tv_usec * (int) pid;
#ifdef MPI
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  seed = seed * np + rank;
#endif
  o->rstr = rnd_alloc(seed);
  if (!o->rstr) {
    free(o->xi);
    free(o);
    return NULL;
  }
  o->lambda = o->xi + n;
  o->upr = o->lambda + n;
  o->tren = o->upr + n;
  o->sig = o->tren + n;
  o->runge_temp = o->sig + n;
  for (i = 0; i < n; ++i) {
    o->xi[i] = xi;
    o->lambda[i] = lambda;
    o->tren[i] = tren;
    o->upr[i] = upr + 2.0 * xi;
  }
  o->upr[0] -= xi;
  o->upr[o->n1] -= xi;
  o->nsite = 0;
  for (i = 0; i < n; ++i) {
    a = o->upr[i] != 0.0;
    o->sig[i] = (real_t) a * _sqrt(o->kt * tren);
    o->nsite += a;
  }

  o->step_autonomic = &runge;
  o->step_coupled = &runge2;
  o->eq_autonomic = &f;
  o->eq_coupled = &f2;
  o->equilibrate = &osc_equilibrate;
  o->x0 = &osc_x0_rand;
  o->del = &osc_free;

  return o;
}

