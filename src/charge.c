#include "charge.h"

#include <stdlib.h>
#include <math.h>

#ifdef USE_MKL
# include <mkl.h>
#else
# include <gsl/gsl_cblas.h>
#endif

extern double a0;
extern double ar[];
extern double ai[];
extern double thr[];
extern double thi[];

#ifdef USE_SINGLE
void thomas1s(int n, const float *ai, float ar, const float *bi, const float *di, const float *dr, float *xi, float *xr, float *work);
float clapack_slamch(char cmach);
void clapack_sstevx(char jobz, char range, int n, float *d, float *e,
                    float vl, float vu, int il, int iu, float abstol,
                    int *m, float *w, float *z, int ldz,
                    float *work, int *iwork, int *ifail, int *info);
#define REAL_MAX FLT_MAX
#else
void thomas1d(int n, const double *ai, double ar, const double *bi, const double *di, const double *dr, double *xi, double *xr, double *work);
double clapack_dlamch(char cmach);
void clapack_dstevx(char jobz, char range, int n, double *d, double *e,
                    double vl, double vu, int il, int iu, double abstol,
                    int *m, double *w, double *z, int ldz,
                    double *work, int *iwork, int *ifail, int *info);
#define REAL_MAX DBL_MAX
#endif

charge_eq_t *mk_charge(int n, int n0, charge_parm_t *def)
{
  charge_eq_t *eq;
  eq = (charge_eq_t *) malloc(sizeof(charge_eq_t) + 11 * n * sizeof(real_t) + 5 * n * sizeof(int));
  if (!eq)
    return NULL;
  eq->d = (real_t *) (eq + 1);
  eq->s = eq->d + n;
  eq->s0 = eq->s + n;
  eq->work = eq->s0 + n;
  eq->n = n;
  eq->n0 = n0;
  eq->chi = def->chi;
  eq->lambda = def->lambda;

  return eq;
}

int charge_init(charge_eq_t *eq, etrans_opt_t *o)
{
  int i;

  if (o->chi->count)
    eq->chi = o->chi->dval[0];
  if (o->lambda->count)
    eq->lambda = o->lambda->dval[0];

  for (i = 0; i < eq->n; i++) {
    eq->d[i] += eq->lambda * (i - eq->n0);
  }
  return 0;
}

void charge_del(charge_eq_t *eq)
{
  free(eq);
}

int charge_write(const charge_eq_t *eq, FILE *f)
{
  int n;

  n = 2 * eq->n - 1;

  if (fwrite(&eq->n0, sizeof(int), 1, f) != 1)
    return -1;
  if (fwrite(&eq->chi, sizeof(real_t), 1, f) != 1)
    return -1;
  if (fwrite(&eq->lambda, sizeof(real_t), 1, f) != 1)
    return -1;
  if (fwrite(eq->d, sizeof(real_t), n, f) != n)
    return -1;

  return 0;
}

int charge_read(charge_eq_t *eq, FILE *f)
{
  int n;
  n = 2 * eq->n - 1;

  if (fread(&eq->n0, sizeof(int), 1, f) != 1)
    return -1;
  if (fread(&eq->chi, sizeof(real_t), 1, f) != 1)
    return -1;
  if (fread(&eq->lambda, sizeof(real_t), 1, f) != 1)
    return -1;
  
  if (fread(eq->d, sizeof(real_t), n, f) != n)
    return -1;

  return 0;
}

real_t charge_expM(charge_eq_t *eq, int m, real_t h0, real_t sn, real_t cs, const real_t *u, const real_t *v, real_t *sr, real_t *si)
{
  real_t *s0r, *s0i, *tr, *ti, *dk, *ds, *work, b2;
  int k, l;
  /* 8n */
  s0r = eq->work;
  s0i = s0r + eq->n;
  tr = s0i + eq->n;
  ti = tr + eq->n;
  ds = ti + eq->n;
  dk = ds + eq->n;
  work = dk + eq->n;

  /* ds = h * (d0 + dt * v(j)) */
  cblas_copy(eq->n, eq->d, 1, ds, 1);
  cblas_axpy(eq->n, eq->chi, u, 1, ds, 1);
  cblas_axpy(eq->n, (m + 0.5) * h0 * eq->chi, v, 1, ds, 1);
  cblas_scal(eq->n, h0, ds, 1);

  cblas_copy(eq->n, sr, 1, s0r, 1);
  cblas_copy(eq->n, si, 1, s0i, 1);

  cblas_scal(eq->n, a0, s0r, 1);
  cblas_scal(eq->n, a0, s0i, 1);


  for (l = 0; l < 14; ++l) {
    for (k = 0; k < eq->n; ++k)
      dk[k] = ds[k] - thi[l];
    thomas1(eq->n, dk, -thr[l], eq->s0, si, sr, ti, tr, work);
    /* s0r = s0r + ar[k] * tr - ai[k] * ti */
    cblas_axpy(eq->n,  ar[l], tr, 1, s0r, 1);
    cblas_axpy(eq->n, -ai[l], ti, 1, s0r, 1);
    /* s0i = s0i + ai[k] * tr + ar[k] * ti */
    cblas_axpy(eq->n,  ai[l], tr, 1, s0i, 1);
    cblas_axpy(eq->n,  ar[l], ti, 1, s0i, 1);
  }

  cblas_copy(eq->n, s0r, 1, sr, 1);
  cblas_copy(eq->n, s0i, 1, si, 1);
  /* sr = cs * s0r + sn * s0i */
  cblas_scal(eq->n, cs, sr, 1);
  cblas_axpy(eq->n, sn, s0i, 1, sr, 1);
  /* si = cs * s0i - sn * s0r */
  cblas_scal(eq->n,  cs, si, 1);
  cblas_axpy(eq->n, -sn, s0r, 1, si, 1);

  b2 = cblas_dot(2 * eq->n, sr, 1, sr, 1);
  /*
  for (k = 0; k < n; ++k) {
    ds[k] = h0 * (eq->d[k] + eq->chi * (u[k] + (m + 0.5) * h0 * v[k]));
    s0r[k] = a0 * sr[k];
    s0i[k] = a0 * si[k];
  }
  for (l = 0; l < 14; ++l) {
    for (k = 0; k < n; ++k)
      dk[k] = ds[k] - thi[l];
    thomas1(eq->n, dk, -thr[l], eq->s0, si, sr, ti, tr, work);
    for (k = 0; k < n; ++k) {
      s0r[k] += ar[l] * tr[k] - ai[l] * ti[k];
      s0i[k] += ai[l] * tr[k] + ar[l] * ti[k];
    }
  }
  b2 = 0.0;
  for (k = 0; k < n; ++k) {
    sr[k] = cs * s0r[k] + sn * s0i[k];
    si[k] = cs * s0i[k] - sn * s0r[k];
    b2 += sr[k] * sr[k] + si[k] * si[k];
  }
  */
  return b2;
}

void charge_step_adj(charge_eq_t *eq, real_t h, int revstep, const real_t *u, const real_t *v, const real_t *sr, const real_t *si,
		     int *nskip, real_t *h0, real_t *sn, real_t *cs, real_t *s0)
{
  static real_t abstol;

  real_t *d0, *ds, *wa, *work;
  real_t w, a, b, hmax;
  int k, info, neigs;
  int *iwork;

  abstol = 2.0 * clapack_lamch('S');

  /* 8n real + 5n int */
  d0 = eq->work;
  ds = d0 + eq->n;
  wa = ds + eq->n;
  work = wa + eq->n;
  iwork = (int *) (work + 5 * eq->n);

  /* d0 = d + u(j) */
  cblas_copy(eq->n, eq->d, 1, d0, 1);
  cblas_axpy(eq->n, eq->chi, u, 1, d0, 1);
  cblas_copy(eq->n, d0, 1, ds, 1);
  cblas_axpy(eq->n, revstep * h * eq->chi, v, 1, ds, 1);

  clapack_stevx('N', 'I', eq->n, d0, eq->s, 0.0, 1.0, 1, 1, abstol, &neigs, wa, NULL, eq->n, work, iwork, NULL, &info);
  a = wa[0];
  clapack_stevx('N', 'I', eq->n, ds, eq->s, 0.0, 1.0, 1, 1, abstol, &neigs, wa, NULL, eq->n, work, iwork, NULL, &info);
  b = wa[0];

  w = a < b ? a : b;
  for (k = 0; k < eq->n; ++k) {
    d0[k] -= w;
    ds[k] -= w;
  }
  clapack_stevx('N', 'I', eq->n, d0, eq->s, 0.0, 1.0, eq->n, eq->n, abstol, &neigs, wa, NULL, eq->n, work, iwork, NULL, &info);
  a = wa[0];
  clapack_stevx('N', 'I', eq->n, ds, eq->s, 0.0, 1.0, eq->n, eq->n, abstol, &neigs, wa, NULL, eq->n, work, iwork, NULL, &info);
  b = wa[0];

  hmax = 0.75 * M_PI / (a > b ? a : b);
  /*if (hmax > d->osc_h) hmax = d->osc_h;*/
  *nskip = h > hmax ? (int) ceil(h / hmax) : 1;
  *h0 = h / *nskip;
  /*
  if (h0 > mx) mx = h0;
  if (h0 < mn) mn = h0;
  */

  a = w * *h0;
  *sn = 0.5 * _sin(a);
  *cs = 0.5 * _cos(a);

  cblas_copy(eq->n - 1, eq->s, 1, s0, 1);
  cblas_scal(eq->n - 1, *h0, s0, 1);
}
