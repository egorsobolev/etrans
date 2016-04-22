#include "sol.h"
#include "chain.h"

#include <math.h>

void int_E(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r)
{
  int n, n1, i;
  double e;
  const real_t *v, *u, *y;

  n = eq->n;
  n1 = n - 1;
  y = x + n;
  u = x + 2 * n;
  v = x + 3 * n;

  e = 0.0;
  for (i = 0; i < n; ++i) {
    e += v[i] * v[i];
  }
  e = 0.5 * e;
  for (i = 0; i < n; ++i) {
    e += (double) eq->ch->chi * u[i] * b2[i];
  }
  for (i = 0; i < n; ++i) {
    e += (double) eq->ch->d[i] * (x[i] * x[i] + y[i] * y[i]);
  }
  for (i = 0; i < n1; ++i) {
    e += (double) 2.0 * eq->ch->s[i] * (x[i] * x[i+1] + y[i] * y[i+1]);
  }

  *r = e + eq->chain->en_potential(eq->chain, u);
}

void int_Eu(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r)
{
  *r = eq->chain->en_potential(eq->chain, x + 2 * eq->n);
}

void int_Ev(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r)
{
  int n, i;
  double e;
  const real_t *v;

  n = eq->n;
  v = x + 3 * n;

  e = 0.0;
  for (i = 0; i < n; ++i) {
    e += v[i] * v[i];
  }
  *r = 0.5 * e;
}

void int_N(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r)
{
  int i;
  double nrm;

  nrm = 0.0;
  for (i = 0; i < eq->n; ++i)
    nrm += b2[i];

  *r = _fabs(1.0 - nrm);
}

void int_S(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r)
{
  int i;
  double s;

  s = 0.0;
  for (i = 0; i < eq->n; ++i) {
    if (b2[i] > 0.0)
      s -= b2[i] * _log(b2[i]);
  }
  *r = s;
}

void int_x(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r)
{
  int i;
  double j, x2;

  x2 = 0.0;
  j = (double) -eq->ch->n0;
  for (i = 0; i < eq->n; ++i) {
    x2 += (double) b2[i] * j;
    j += 1.0;
  }
  *r = x2;
}

void int_dx(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r)
{
  int n, i;
  double j, x2;

  const real_t *y, *dy;

  n = eq->n;
  y = x + n;
  dy = dx + n;

  x2 = 0.0;
  j = (double) -eq->ch->n0;
  for (i = 0; i < n; ++i) {
    x2 += (double) 2.0 * (x[i] * dx[i] + y[i] * dy[i]) * j;
    j += 1.0;
  }
  *r = x2;
}

void int_x2(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r)
{
  int i;
  double j, x2;

  x2 = 0.0;
  j = (double) -eq->ch->n0;
  for (i = 0; i < eq->n; ++i) {
    x2 += (double) b2[i] * j * j;
    j += 1.0;
  }
  *r = x2;
}

void int_dx2(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r)
{
  int n, i;
  double j, x2;
  const real_t *y, *dy;

  n = eq->n;
  y = x + n;
  dy = dx + n;

  x2 = 0.0;
  j = (double) -eq->ch->n0;
  for (i = 0; i < n; ++i) {
    x2 += (double) 2.0 * (x[i] * dx[i] + y[i] * dy[i]) * j * j;
    j += 1.0;
  }
  *r = x2;
}

void int_d2x2(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r)
{
  int n, i;
  double j, x2;
  const real_t *y, *dy, *d2y;

  n = eq->n;
  y = x + n;
  dy = dx + n;
  d2y = d2x + n;

  x2 = 0.0;
  j = (double) -eq->ch->n0;
  for (i = 0; i < n; ++i) {
    x2 += (double) 2.0 * (dx[i] * dx[i] + dy[i] * dy[i] + x[i] * d2x[i] + y[i] * d2y[i]) * j * j;
    j += 1.0;
  }
  *r = x2;
}

void int_j(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r)
{
  int n, n1, i;
  const real_t *y;
  double e;
  n = eq->n;
  n1 = n - 1;
  y = x + n;
  e = 0.0;
  for (i = 0; i < n1; ++i) {
    e += (double) x[i] * y[i+1] - y[i] * x[i+1];
  }
  *r = -2.0 * e;
}

void int_Eq(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r)
{
  int n, n1, i;
  const real_t *y;
  double e;
  n = eq->n;
  n1 = n - 1;
  y = x + n;
  e = 0.0;
  for (i = 0; i < n; ++i) {
    e += (double) eq->ch->d[i] * (x[i] * x[i] + y[i] * y[i]);
  }
  for (i = 0; i < n1; ++i) {
    e += (double) 2.0 * eq->ch->s[i] * (x[i] * x[i+1] + y[i] * y[i+1]);
  }
  *r = e;
}

void int_Eb(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r)
{
  int i, n;
  const real_t *u;
  double e;
  n = eq->n;
  u = x + 2 * n;
  e = 0.0;
  for (i = 0; i < n; ++i) {
    e += (double) eq->ch->chi * u[i] * b2[i];
  }
  *r = e;
}

void int_invb4(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r)
{
  int i, n;
  double e;
  n = eq->n;
  e = 0.0;
  for (i = 0; i < n; ++i) {
    e += (double) b2[i] * b2[i];
  }
  *r = 1.0 / e;
}

void cp_u(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r)
{
  int i;
  const real_t *u;

  u = x + 2 * eq->n;
  for (i = 0; i < eq->n; i++)
    r[i] = (double) u[i];
}

void cp_v(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r)
{
  int i;
  const real_t *v;

  v = x + 3 * eq->n;
  for (i = 0; i < eq->n; i++)
    r[i] = (double) v[i];
}

void cp_b2(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r)
{
  int i;
  for (i = 0; i < eq->n; i++)
    r[i] = (double) b2[i];
}



int len_chain(int n)
{
  return n;
}
int len_scalar(int n)
{
  return 1;
}

long size_chain_traj(int n, int nstep)
{
  return n * nstep;
}
long size_scalar_traj(int n, int nstep)
{
  return nstep;
}
long size_chain_scalar(int n, int nstep)
{
  return n;
}

int_func_t int_func[MAX_FUNC] = {
  /* 01 */ "x", &int_x, &size_scalar_traj, &len_scalar,
  /* 02 */ "dx", &int_dx, &size_scalar_traj, &len_scalar,
  /* 03 */ "x2", &int_x2, &size_scalar_traj, &len_scalar,
  /* 04 */ "dx2", &int_dx2, &size_scalar_traj, &len_scalar,
  /* 05 */ "d2x2", &int_d2x2, &size_scalar_traj, &len_scalar,
  /* 06 */ "Eq", &int_Eq, &size_scalar_traj, &len_scalar,
  /* 07 */ "Eb", &int_Eb, &size_scalar_traj, &len_scalar,
  /* 08 */ "Eu", &int_Eu, &size_scalar_traj, &len_scalar,
  /* 09 */ "Ev", &int_Ev, &size_scalar_traj, &len_scalar,
  /* 10 */ "E", &int_E, &size_scalar_traj, &len_scalar,
  /* 11 */ "invb4", &int_invb4, &size_scalar_traj, &len_scalar,
  /* 12 */ "j", &int_j, &size_scalar_traj, &len_scalar,
  /* 13 */ "S", &int_S, &size_scalar_traj, &len_scalar,
  /* 14 */ "N", &int_N, &size_scalar_traj, &len_scalar,
  /* 15 */ "p", &cp_b2, &size_chain_traj, &len_chain,
  /* 16 */ "u", &cp_u, &size_chain_traj, &len_chain,
  /* 17 */ "v", &cp_v, &size_chain_traj, &len_chain,
};
