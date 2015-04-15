#include "config.h"
#include "osc.h"

#include <time.h>
#include <string.h>
#include <math.h>

charge_parm_t osc_charge_defs = {
  0.02, 0.0, /* chi_e */
  0.0, 6.84, 10.0, 10.0, /*  G  A  T  C */
  1.276, 1.352, 2.081, 1.671, /* GG GA GT GC */
  0.744, 0.456, 1.595, 0.927, /* AG AA AT AC */
  1.291, 1.307, 2.400, 1.155, /* TG TA TT TC */
  0.638, 0.441, 1.519, 0.623  /* CG CA CT CC */
};

real_t en_potential(const chain_eq_t *chain, const real_t *u)
{
  int k;
  real_t e0, e1;
  const oscilator_t *o = (oscilator_t *) chain;
  e0 = 0.0;
  for (k = 0; k < o->n; ++k) {
    e0 += u[k] * u[k];
  }
  e1 = 0.0;
  for (k = 1; k < o->n; ++k) {
    e1 += u[k-1] * u[k];
  }
  return 0.5 * (o->K + o->D) * e0 - o->D * e1;
}

void osc_dv(const chain_eq_t *chain, const real_t *x, real_t *dv)
{
  int i;
  const oscilator_t *o = (oscilator_t *) chain;
  const real_t *v;
  v = x + o->n;

  dv[0] = -(o->K + o->D) * x[0] - o->L * v[0] + o->D * x[1];
  for (i = 1; i < o->n1; ++i)
    dv[i] = -(o->K + 2.0 * o->D) * x[i] - o->L * v[i] + o->D * (x[i + 1] + x[i - 1]);
  dv[o->n1] = -(o->K + o->D) * x[o->n1] - o->L * v[o->n1] + o->D * x[o->n1 - 1];
}

void osc_holstein_dv(const chain_eq_t *chain, const real_t *x, real_t *dv, real_t *b2)
{
  int i;
  const oscilator_t *o = (oscilator_t *) chain;
  const real_t *v;
  v = x + o->n;

  dv[0] = -(o->K + o->D) * x[0] - o->L * v[0] + o->D * x[1] - o->mu * b2[0];
  for (i = 1; i < o->n1; ++i)
    dv[i] = -(o->K + 2.0 * o->D) * x[i] - o->L * v[i] + o->D * (x[i + 1] + x[i - 1]) - o->mu * b2[i];
  dv[o->n1] = -(o->K + o->D) * x[o->n1] - o->L * v[o->n1] + o->D * x[o->n1 - 1] - o->mu * b2[o->n1];
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
  fcnt = fwrite(&c->sigF, sizeof(real_t), 6, f);
  if (fcnt != 6)
    return -1;

  return 0;
}

int osc_read(chain_eq_t *chain, FILE *f)
{
  oscilator_t *c = (oscilator_t *) chain;
  size_t fcnt;

  fcnt = fread(&c->sigF, sizeof(real_t), 6, f);
  if (fcnt != 6)
    return -1;

  return 0;
}

size_t osc_nbytes(int n)
{
  return sizeof(oscilator_t);
}

oscilator_t *mk_osc(int n)
{
  static real_t e0 = 0.261838952;
  static real_t t0 = 100.0;
  oscilator_t *c;

  c = (oscilator_t *) malloc(osc_nbytes(n));
  if (!c)
    return NULL;

  c->n = n;
  c->n1 = n - 1;

  c->dv = &osc_dv;
  c->dv_hst = &osc_holstein_dv;
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
