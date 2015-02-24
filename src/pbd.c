#include "config.h"
#include "pbd.h"

#include <time.h>
#include <string.h>
#include <math.h>

#ifdef USE_MKL
# include <mkl.h>
#else
# include <gsl/gsl_cblas.h>
#endif

charge_parm_t pbd_charge_defs = {
  8.809624572754315, 0.0,/* chi_e */
  264.4117479078003, 360.3676241646633, 405.1470330845326, 405.1470330845326, /*  G  A  T  C */
   17.9117635679478,  18.9779399708018,  29.2132334382005,  23.4558808627887, /* GG GA GT GC */
   10.4485287479695,   6.3970584171242,  22.3897044599347,  13.0073521148192, /* AG AA AT AC */
   18.1249988485186,  18.3382341290894,  33.6911743301875,  16.2058813233813, /* TG TA TT TC */
    8.9558817839739,   6.1838231365534,  21.3235280570807,   8.7426465034031  /* CG CA CT CC */
};

real_t ekin_pbd(int n, const real_t *v)
{
  real_t e;
  e = 0.5 * cblas_dot(n, v, 1, v, 1);
  return e;
}

real_t epot_pbd(const chain_eq_t *chain, const real_t *u)
{
  real_t e1, e2, a;
  int i;
  const pbd_chain_t *c = (pbd_chain_t *) chain;
  e1 = 0.0;
  for (i = 0; i < c->n; i++) {
    a = exp(-c->sigma * u[i]) - 1.0;
    e1 += a * a;
  }
  e2 = 0.0;
  for (i = 0; i < c->n1; i++) {
    a = u[i+1] - u[i];
    e2 += (1.0 + c->rho * exp(-c->epsilon * (u[i] + u[i+1]))) * a * a;
  }
  return 0.5 * (c->omegaM2 * e1 + c->omegaB2 * e2);
}

void eq_pbd(const chain_eq_t *chain, const real_t *x, real_t *dx)
{
  real_t *dv;
  const real_t *v;
  int i;
  real_t A, Er, El, Dr, Dl;
  const pbd_chain_t *c = (pbd_chain_t *) chain;

  dv = dx + c->n;
  v = x + c->n;
  cblas_copy(c->n, v, 1, dx, 1);

  A = exp(-c->sigma * x[0]);
  Dr = x[0] - x[1];
  Er = c->rho * exp(-c->epsilon * (x[0] + x[1])); 
  dv[0] = c->omegaM2 * (A - 1.0) * A +
    c->omegaB2 * Dr * (0.5 * c->epsilon * (Dr * Er) - Er - 1.0)
    - c->Gamma * v[0];

  for (i = 1; i < c->n1; ++i) {
    A = exp(-c->sigma * x[i]);
    Dl = x[i] - x[i-1];
    Dr = x[i] - x[i+1];
    El = c->rho * exp(-c->epsilon * (x[i-1] + x[i]));
    Er = c->rho * exp(-c->epsilon * (x[i+1] + x[i])); 
    dv[i] = c->omegaM2 * (A - 1.0) * A + c->omegaB2 * 
      (
       Dl * (0.5 * c->epsilon * (Dl * El) - El - 1.0) +
       Dr * (0.5 * c->epsilon * (Dr * Er) - Er - 1.0)
       ) - c->Gamma * v[i];
  }
  A = exp(-c->sigma * x[c->n1]);
  Dl = x[c->n1] - x[c->n1-1];
  El = c->rho * exp(-c->epsilon * (x[c->n1-1] + x[c->n1])); 
  dv[c->n1] = c->omegaM2 * (A - 1.0) * A +
    c->omegaB2 * Dl * (0.5 * c->epsilon * (Dl * El) - El - 1.0)
    - c->Gamma * v[c->n1];
}

void eq_pbd_holstein(const chain_eq_t *chain, const real_t *x, real_t *dx, real_t *b2)
{
  real_t *dv;
  const real_t *v;
  int i;
  real_t A, Er, El, Dr, Dl;
  const pbd_chain_t *c = (pbd_chain_t *) chain;

  dv = dx + c->n;
  v = x + c->n;
  cblas_copy(c->n, v, 1, dx, 1);

  A = exp(-c->sigma * x[0]);
  Dr = x[0] - x[1];
  Er = c->rho * exp(-c->epsilon * (x[0] + x[1])); 
  dv[0] = c->omegaM2 * (A - 1.0) * A + 
    c->omegaB2 * Dr * (0.5 * c->epsilon * (Dr * Er) - Er - 1.0)
    - c->Gamma * v[0] - c->chi * b2[0];

  for (i = 1; i < c->n1; ++i) {
    A = exp(-c->sigma * x[i]);
    Dl = x[i] - x[i-1];
    Dr = x[i] - x[i+1];
    El = c->rho * exp(-c->epsilon * (x[i-1] + x[i]));
    Er = c->rho * exp(-c->epsilon * (x[i+1] + x[i])); 
    dv[i] = c->omegaM2 * (A - 1.0) * A + c->omegaB2 * 
      (Dl * (0.5 * c->epsilon * (Dl * El) - El - 1.0) +
       Dr * (0.5 * c->epsilon * (Dr * Er) - Er - 1.0))
      - c->Gamma * v[i] - c->chi * b2[i];
  }
  A = exp(-c->sigma * x[c->n1]);
  Dl = x[c->n1] - x[c->n1-1];
  El = c->rho * exp(-c->epsilon * (x[c->n1-1] + x[c->n1])); 
  dv[c->n1] = c->omegaM2 * (A - 1.0) * A + 
    c->omegaB2 * Dl * (0.5 * c->epsilon * (Dl * El) - El - 1.0)
    - c->Gamma * v[c->n1] - c->chi * b2[c->n1];
}

void runge_pbd(chain_eq_t *chain, real_t h, real_t *x)
{
  real_t sqrth, h_2, *y, *k1, *k2;
  int i, m;
  pbd_chain_t *c = (pbd_chain_t *) chain;

  sqrth = _sqrt(h) * c->sigmaF;
  h_2 = 0.5 * h;
  m = 2 * c->n;
  y = c->runge_temp;
  k1 = y + m;
  k2 = k1 + m;

  cblas_copy(m, x, 1, y, 1);

  rnd_gaussian(c->n, x + c->n, 1.0);
  /* memset(x + o->n, 0, o->n * sizeof(float)); */

  for (i = c->n; i < m; ++i)
    x[i] = sqrth * x[i] + y[i];
  eq_pbd(chain, x, k1);
  cblas_axpy(m, h, k1, 1, y, 1);
  eq_pbd(chain, y, k2);
  cblas_axpy(m, h_2, k1, 1, x, 1);
  cblas_axpy(m, h_2, k2, 1, x, 1);
}

void runge_pbd_holstein(chain_eq_t *chain, real_t h, real_t *x, real_t *b0, real_t *b1)
{
  real_t sqrth, h_2, *y, *k1, *k2;
  int i, m;
  pbd_chain_t *c = (pbd_chain_t *) chain;

  sqrth = _sqrt(h) * c->sigmaF;
  h_2 = 0.5 * h;
  m = 2 * c->n;
  y = c->runge_temp;
  k1 = y + m;
  k2 = k1 + m;

  cblas_copy(m, x, 1, y, 1);

  rnd_gaussian(c->n, x + c->n, 1.0);
  /*memset(x + o->n, 0, o->n * sizeof(float));*/

  for (i = c->n; i < m; ++i)
    x[i] = sqrth * x[i] + y[i];
  eq_pbd_holstein(chain, x, k1, b0);
  cblas_axpy(m, h, k1, 1, y, 1);
  eq_pbd_holstein(chain, y, k2, b1);
  cblas_axpy(m, h_2, k1, 1, x, 1);
  cblas_axpy(m, h_2, k2, 1, x, 1);
}

void pbd_equilibrate(chain_eq_t *chain)
{
  pbd_chain_t *c = (pbd_chain_t *) chain;
  int i, j, n, m, s0, nskip, nr, l, k;
  real_t *u, *v, ep, ek, ep2, ek2, a, b;
  FILE *f;
  real_t *t, *r, *Mu, *Du;
  long pos;

  u = c->f0;
  v = u + c->n;

  s0 = c->n / 2;
  nr = 500;

  Mu = (real_t *) malloc(2 * (nr + 1) * c->n * sizeof(real_t));
  Du = Mu + c->n;
  r = Du + c->n;
  t = r + nr * c->n;

  f = fopen("traj.dat", "wb");
  fwrite(&c->n, sizeof(int), 1, f);
  pos = ftell(f);
  fwrite(&m, sizeof(int), 1, f);

  /* heating */
  ek2 = 0.5 * cblas_dot(c->n, v, 1, v, 1) / c->n;

  m = 0;
  i = 0;
  while (i < 1000) {
    runge_pbd(chain, c->h, c->f0);
    ek = 0.5 * cblas_dot(c->n, v, 1, v, 1) / c->n;
    i += ((2.0 * ek2 - c->theta) * (2.0 * ek - c->theta)) < 0.0;

    ep = epot_pbd(chain, u) / c->n;

    fwrite(&ek, sizeof(real_t), 1, f);
    fwrite(&ep, sizeof(real_t), 1, f);

    ek2 = ek;
    m++;
  };

  nskip = 200;
  /* autocorrelation function */
  for (i = nr - 1; i >= 0; i--) {
    for (j = 0; j < nskip; j++)
      runge_pbd(chain, c->h, c->f0);
    memcpy(t + i * c->n, u, c->n * sizeof(real_t)); 
    /*t[i] = u[s0]*/;
  }
  /*Mu = Du = 0.0;*/
  memset(Mu, 0, c->n * sizeof(real_t));
  memset(Du, 0, c->n * sizeof(real_t));
  memset(r, 0, nr * c->n * sizeof(real_t));
  for (i = 0; i < 10000; i++) {
    for (j = 0; j < nskip; j++)
      runge_pbd(chain, c->h, c->f0);
    for (j = 0; j < c->n; j++) {
      Mu[j] += u[j];
      Du[j] += u[j] * u[j];
    }
    l = nr * c->n - 1;
    for (j = nr - 1; j > 0; j--) {
      for (k = c->n - 1; k >= 0; k--) {
	r[l] += u[k] * t[l];
	t[l] = t[l - c->n];
	l--;
      }
      /*
	r[j] += cblas_dot(c->n, u, 1, t + j * c->n, 1);
	memcpy(t + j * c->n, t + j * c->n - c->n, c->n * sizeof(real_t));
      */
    }
    for (k = c->n - 1; k >=0; k--) {
      r[k] += u[k] * t[k];
      t[k] = u[k];
    }
    /*
    r[0] += cblas_dot(c->n, u, 1, t, 1);
    memcpy(t, u, sizeof(real_t));
    */

    /*
    a = u[s0];
    Mu += a;
    Du += a * a;
    for (j = 99; j > 0; j--) {
      r[j] += a * t[j];
      t[j] = t[j-1];
    }
    r[0] += a * t[0];
    t[0] = a;
    */
  }

  for (k = 0; k < c->n; k++) {
    Mu[k] = Mu[k] * Mu[k] / 10000.0;
    Du[k] -= Mu[k];
  }
  l = 0;
  for (j = 0; j < nr; j++)
    for (k = 0; k < c->n; k++) {
      r[l] = (r[l] - Mu[k]) / Du[k];
      l++;
    }

  /*
  a = Mu * Mu / (10000.0 * c->n);
  Du -= a;
  for (j = 0; j < nr; j++)
    r[j] = (r[j] - a) / Du;
  */

  fwrite(r, sizeof(real_t), nr * c->n, f);

  fseek(f, pos, SEEK_SET);
  fwrite(&m, sizeof(int), 1, f);

  fclose(f);

  free(Mu);
}

void pbd_x0_rand(chain_eq_t *chain, real_t *x)
{
  pbd_chain_t *c = (pbd_chain_t *) chain;
  int i, n;

  n = 50000;
  for (i = 0; i < n; ++i)
    runge_pbd(chain, c->h, c->f0);

  memcpy(x, c->f0, 2 * c->n * sizeof(real_t));
}

void pbd_free(chain_eq_t *chain)
{
  free(chain);
}

int pbd_init(chain_eq_t *chain, etrans_opt_t *o)
{
  pbd_chain_t *c = (pbd_chain_t *) chain;

  if (o->gamma->count)
    c->Gamma = o->gamma->dval[0];
  c->theta0 = 1.07716655e-3;
  if (o->temp->count)
    c->theta = c->theta0 * o->temp->dval[0];
  c->sigmaF = _sqrt(2.0 * c->theta * c->Gamma);

  if (o->omegaM2->count)
    c->omegaM2 = o->omegaM2->dval[0]; /* 1.0 */
  if (o->sigma->count)
    c->sigma = o->sigma->dval[0]; /* 1.0 */

  if (o->omegaB2->count)
    c->omegaB2 = o->omegaB2->dval[0]; /* 2.5249337204898369e-2 */
  if (o->rho->count)
    c->rho = o->rho->dval[0]; /* 0.5 */
  if (o->epsilon->count)
    c->epsilon = o->epsilon->dval[0]; /* 7.8651685393258425e-2 */

  if (o->mu->count)
    c->chi = o->mu->dval[0];

  if (o->h->count)
    c->h = o->h->dval[0];
  else
    c->h = 0.005 / _sqrt(c->omegaM2 > c->omegaB2 ? c->omegaM2 : c->omegaB2); 

  return 0;
}

int pbd_write(const chain_eq_t *chain, FILE *f)
{
  const pbd_chain_t *c = (pbd_chain_t *) chain;
  size_t fcnt;

  fcnt = fwrite("PBD", sizeof(char), 3, f);
  if (fcnt != 3)
    return -1;
  fcnt = fwrite(&c->omegaM2, sizeof(real_t), 10, f);
  if (fcnt != 10)
    return -1;

  return 0;
}

int pbd_read(chain_eq_t *chain, FILE *f)
{
  pbd_chain_t *c = (pbd_chain_t *) chain;
  size_t fcnt;

  fcnt = fread(&c->omegaM2, sizeof(real_t), 10, f);
  if (fcnt != 10)
    return -1;

  return 0;
}

size_t pbd_nbytes(int n)
{
  return sizeof(pbd_chain_t) + 8 * n * sizeof(real_t);
}

pbd_chain_t *mk_pbd(int n)
{
  pbd_chain_t *c;
  int m;

  /* 6n - runge workspace */
  /* 2n - stored termolized x0 */
  c = (pbd_chain_t *) malloc(pbd_nbytes(n));
  if (!c)
    return NULL;

  c->f0 = (real_t *) (c + 1);
  c->runge_temp = c->f0 + m;

  c->n = n;
  c->n1 = n - 1;
  m = 2 * c->n;

  memset(c->f0, 0, m * sizeof(real_t));

  c->Gamma = 0.084212403085279;
  c->theta0 = 1.07716655e-3;
  c->theta = c->theta0 * 300;
  c->sigmaF = _sqrt(2.0 * c->theta * c->Gamma);

  c->omegaM2 = 1.0;
  c->sigma = 1.0;

  c->omegaB2 = 2.5249337204898369e-2;
  c->rho = 0.5;
  c->epsilon = 7.8651685393258425e-2;

  c->chi = 0.516426300866580;

  c->h = 0.005 / _sqrt(c->omegaM2 > c->omegaB2 ? c->omegaM2 : c->omegaB2); 


  c->step_autonomic = &runge_pbd;
  c->step_coupled = &runge_pbd_holstein;
  c->eq_autonomic = &eq_pbd;
  c->eq_coupled = &eq_pbd_holstein;
  c->equilibrate = &pbd_equilibrate;
  c->x0 = &pbd_x0_rand;
  c->del = &pbd_free;
  c->en_potential = &epot_pbd;
  c->init = &pbd_init;
  c->write = &pbd_write;
  c->read = &pbd_read;
  c->nbytes = &pbd_nbytes;

  return c;
}
