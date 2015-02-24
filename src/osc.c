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

charge_parm_t osc_charge_defs = {
  0.02, 0.0, /* chi_e */
  0.0, 6.84, 10.0, 10.0, /*  G  A  T  C */
  1.276, 1.352, 2.081, 1.671, /* GG GA GT GC */
  0.744, 0.456, 1.595, 0.927, /* AG AA AT AC */
  1.291, 1.307, 2.400, 1.155, /* TG TA TT TC */
  0.638, 0.441, 1.519, 0.623  /* CG CA CT CC */
};

real_t en_kinetic(int n, const real_t *v)
{
  real_t e;
  e = 0.5 * cblas_dot(n, v, 1, v, 1);
  return e;
}

real_t en_potential(const chain_eq_t *chain, const real_t *u)
{
  real_t e;
  const oscilator_t *o = (oscilator_t *) chain;
  e = 0.5 * (o->K + o->D) * cblas_dot(o->n, u, 1, u, 1) - o->D * cblas_dot(o->n1, u, 1, u + 1, 1);
  return e;
}

void f(const chain_eq_t *chain, const real_t *x, real_t *dx)
{
  real_t *dv;
  const real_t *v;
  int i;
  const oscilator_t *o = (oscilator_t *) chain;

  dv = dx + o->n;
  v = x + o->n;
  cblas_copy(o->n, v, 1, dx, 1);

  dv[0] = -(o->K + o->D) * x[0] - o->L * v[0] + o->D * x[1];
  for (i = 1; i < o->n1; ++i)
    dv[i] = -(o->K + 2.0 * o->D) * x[i] - o->L * v[i] + o->D * (x[i + 1] + x[i - 1]);
  dv[o->n1] = -(o->K + o->D) * x[o->n1] - o->L * v[o->n1] + o->D * x[o->n1 - 1];
  /*
  dv[0] = -o->upr[0] * x[0] - o->tren[0] * v[0] + o->xi[0] * x[1];
  for (i = 1; i < o->n1; ++i)
    dv[i] = -o->upr[i] * x[i] - o->tren[i] * v[i] + o->xi[i] * (x[i + 1] + x[i - 1]);
  dv[o->n1] = -o->upr[o->n1] * x[o->n1] - o->tren[o->n1] * v[o->n1] + o->xi[o->n1] * x[o->n1 - 1];
  */
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
  dv[0] = -(o->K + o->D) * x[0] - o->L * v[0] + o->D * x[1] - o->mu * b2[0];
  for (i = 1; i < o->n1; ++i)
    dv[i] = -(o->K + 2.0 * o->D) * x[i] - o->L * v[i] + o->D * (x[i + 1] + x[i - 1]) - o->mu * b2[i];
  dv[o->n1] = -(o->K + o->D) * x[o->n1] - o->L * v[o->n1] + o->D * x[o->n1 - 1] - o->mu * b2[o->n1];
  /*
  dv[0] = -o->upr[0] * x[0] - o->tren[0] * v[0] + o->xi[0] * x[1] - o->lambda[0] * b2[0];
  for (i = 1; i < o->n1; ++i)
    dv[i] = -o->upr[i] * x[i] - o->tren[i] * v[i] + o->xi[i] * (x[i + 1] + x[i - 1]) - o->lambda[i] * b2[i];
  dv[o->n1] = -o->upr[o->n1] * x[o->n1] - o->tren[o->n1] * v[o->n1] + o->xi[o->n1] * x[o->n1 - 1] - o->lambda[o->n1] * b2[o->n1];
  */
}

void runge(chain_eq_t *chain, real_t h, real_t *x)
{
  real_t sqrth, h_2, *y, *k1, *k2;
  int i, m;
  oscilator_t *o = (oscilator_t *) chain;

  /*
  sqrth = _sqrt(h);
  */
  sqrth = _sqrt(h) * o->sigF;
  h_2 = 0.5 * h;
  m = 2 * o->n;
  y = o->runge_temp;
  k1 = y + m;
  k2 = k1 + m;

  cblas_copy(m, x, 1, y, 1);

  rnd_gaussian(o->n, x + o->n, 1.0);
  /* memset(x + o->n, 0, o->n * sizeof(float)); */

  for (i = o->n; i < m; ++i)
    x[i] = sqrth * x[i] + y[i];
  /*
    x[i] = sqrth * o->sig[i - o->n] * x[i] + y[i];
  */
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

  /*
  sqrth = _sqrt(h);
  */
  sqrth = _sqrt(h) * o->sigF;
  h_2 = 0.5 * h;
  m = 2 * o->n;
  y = o->runge_temp;
  k1 = y + m;
  k2 = k1 + m;

  cblas_copy(m, x, 1, y, 1);

  rnd_gaussian(o->n, x + o->n, 1.0);
  /*memset(x + o->n, 0, o->n * sizeof(float));*/

  for (i = o->n; i < m; ++i)
    x[i] = sqrth * x[i] + y[i];
  /*
    x[i] = sqrth * o->sig[i - o->n] * x[i] + y[i];
  */
  f2(chain, x, k1, b0);
  cblas_axpy(m, h, k1, 1, y, 1);
  f2(chain, y, k2, b1);
  cblas_axpy(m, h_2, k1, 1, x, 1);
  cblas_axpy(m, h_2, k2, 1, x, 1);
}
void osc_equilibrate(chain_eq_t *chain)
{
}

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
  o->sv = _sqrt(o->kt);
  rnd_gaussian(o->n, x + o->n, o->sv);

  if (o->D == 0.0) {
    o->su = _sqrt(o->kt / o->K);
    b = 1.0;
    a = 0.0;
  } else {
    b = 2.0 * o->D;
    a = _sqrt(o->K * (o->K + 2.0 * b)); 

    o->su = _sqrt(o->kt / a);
    xi = (o->K + b - a) / b;

    a = _sqrt(xi);
    b = _sqrt(1.0 - xi);
  }
  rnd_gaussian(o->n, x, o->su);

  if (o->D == 0.0)
    return;

  for (i = 1; i < o->n; ++i) {
    x[i] *= b;
    x[i] += x[i-1] * a; 
  }
}
void osc_free(chain_eq_t *chain)
{
  free(chain);
}

int osc_init(chain_eq_t *chain, etrans_opt_t *o)
{
  static real_t e0 = 0.261838952;
  static real_t t0 = 100.0;
  oscilator_t *c = (oscilator_t *) chain;
  if (o->omegaM2->count)
    c->K = o->omegaM2->dval[0];
  if (o->omegaB2->count)
    c->D = o->omegaB2->dval[0];
  if (o->gamma->count)
    c->L = o->gamma->dval[0];
  if (o->mu->count)
    c->mu = o->mu->dval[0];
  else
    c->mu = o->chi->dval[0];
  if (o->temp->count)
    c->kt = 0.5 * e0 / t0 * o->temp->dval[0];
  c->sigF = _sqrt(2.0 * c->L * c->kt);

  if (o->h->count)
    c->h = o->h->dval[0];
  else
    c->h = 0.01 / _sqrt(c->K > c->D ? c->K : c->D);

  return 0;
}

int osc_write(const chain_eq_t *chain, FILE *f)
{
  oscilator_t *c = (oscilator_t *) chain;
  size_t fcnt;

  fcnt = fwrite("OSC", sizeof(char), 3, f);
  if (fcnt != 3)
    return -1;
  fcnt = fwrite(&c->kt, sizeof(real_t), 6, f);
  if (fcnt != 6)
    return -1;

  return 0;
}

int osc_read(chain_eq_t *chain, FILE *f)
{
  oscilator_t *c = (oscilator_t *) chain;
  size_t fcnt;

  fcnt = fread(&c->kt, sizeof(real_t), 6, f);
  if (fcnt != 6)
    return -1;

  return 0;
}

size_t osc_nbytes(int n)
{
  return sizeof(oscilator_t) + 6 * n * sizeof(real_t);
}

oscilator_t *mk_osc(int n)
{
  static real_t e0 = 0.261838952;
  static real_t t0 = 100.0;
  oscilator_t *c;

  c = (oscilator_t *) malloc(osc_nbytes(n));
  if (!c)
    return NULL;

  c->runge_temp = (real_t *) (c + 1);

  c->n = n;
  c->n1 = n - 1;

  c->step_autonomic = &runge;
  c->step_coupled = &runge2;
  c->eq_autonomic = &f;
  c->eq_coupled = &f2;
  c->equilibrate = &osc_equilibrate;
  c->x0 = &osc_x0_rand;
  c->del = &osc_free;
  c->en_potential = &en_potential;
  c->init = &osc_init;
  c->write = &osc_write;
  c->read = &osc_read;
  c->nbytes = &osc_nbytes;

  c->K = 1e-4;
  c->D = 0.0;
  c->L = 6e-3;
  c->mu = 0.02;
  c->kt = 0.5 * e0 / t0 * 300.0;
  c->sigF = _sqrt(2.0 * c->L * c->kt);

  c->h = 0.01 / _sqrt(c->K > c->D ? c->K : c->D);

  return c;
}
