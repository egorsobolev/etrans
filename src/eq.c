#include "config.h"
#include "etrans.h"
#include "osc.h"
#include "sol.h"

#include <memory.h>
#include <malloc.h>
#include <stdio.h>

#include <math.h>
#include <float.h>

#ifdef USE_MKL
# include <mkl.h>
#else
# include <gsl/gsl_cblas.h>
#endif

#ifdef MPI
#include <mpi.h>
#define walltime() MPI_Wtime()
#else
#include <time.h>
#define walltime()	((double) clock() / CLOCKS_PER_SEC)
#endif

#ifdef USE_SINGLE
#define REAL_MAX FLT_MAX
#else
#define REAL_MAX DBL_MAX
#endif


int eq(eqdata_t *d, real_t h, real_t tmax, int nsamp, sol_t *res)
{
  int nstep, nskip, n, i, j, ln, k, m, q, qn;
  real_t p2, a, sn, cs, h0, normp, mn, mx, nr;
  real_t *u, *v, *b2l, *b2r, *b2t;
  real_t *sr, *si, *dsr, *dsi, *d2sr, *d2si;
  real_t env, enu;
  int rank;
  double elapse, wtm0, wtm1;

#ifdef MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  rank = 0;
#endif


  qn = nsamp;
  n = d->n;
  nstep = d->nstep;

  sr = (real_t *) malloc(10 * n * sizeof(real_t));
  if (!sr) return -1;

  si = sr + n;
  u = si + n;
  v = u + n;
  dsr = v + n;
  dsi = dsr + n;
  d2sr = dsi + n;
  d2si = d2sr + n;
  b2l = d2si + n;
  b2r = b2l + n;

  sol_setzero(res);
   
  elapse = 0;
  wtm0 = walltime();
  for (q = 0; q < qn; q++) {

    d->chain->x0(d->chain, u);

    if (d->x0) {
      cblas_axpy(2 * n, 1.0, d->x0 + 2 * n, 1, u, 1); 
      cblas_copy(2 * n, d->x0, 1, sr, 1);
    } else {
      memset(sr, 0, 2 * n * sizeof(real_t));
      rnd_uniform(1, &a, 2.0 * M_PI);
      sr[d->ch->n0] = cos(a);
      si[d->ch->n0] = sin(a);
    }
    p2 = 1.0;

    for (k = 0; k < n; ++k) {
      b2l[k] = sr[k] * sr[k] + si[k] * si[k];
    }
    diff_b(d->ch, sr, dsr, d2sr);
    sol_reset(res);
    sol_update(d, sr, b2l, dsr, d2sr, res);

    j = d->nskip - 1;
    nr = 0.0;
    mn = REAL_MAX;
    mx = 0.0;
    ln = 0;

    /* cycle by classical steps */
    i = 0;
    while (i < nstep) {

      if (!ln) {
	charge_step_adj(d->ch, h, d->h_revstep, u, v, sr, si, &nskip, &h0, &sn, &cs, d->ch->s0);
	if (h0 > mx) mx = h0;
	if (h0 < mn) mn = h0;
			   
	ln = d->h_revstep;
      }
      for (m = 0; m < nskip; ++m) {
	p2 = charge_expM(d->ch, m, h0, sn, cs, u, v, sr, si);
	normp = _fabs(p2 - 1.0);
	if (normp > nr)
	  nr = normp;

	p2 = _sqrt(p2);
	a = 1.0 / p2;
	cblas_scal(2 * n, a, sr, 1);

      }
      for (k = 0; k < n; ++k) {
	b2r[k] = sr[k] * sr[k] + si[k] * si[k];
      }
      d->chain->step_coupled(d->chain, h, u, b2l, b2r);
      b2t = b2l;
      b2l = b2r;
      b2r = b2t;

      if (!j) {
	diff_b(d->ch, sr, dsr, d2sr);
	sol_update(d, sr, b2l, dsr, d2sr, res);
	j = d->nskip;
	++i;
      }
      --j;
      --ln;
    }
    wtm1 = walltime();
    if (!rank) {
      env = 0.5 * cblas_dot(n, v, 1, v, 1);
      enu = d->chain->en_potential(d->chain, u);
      printf("%6d %8.2f %9.2g %9.6g %9.6g %9.6g %9.6g\n", q + 1, wtm1 - wtm0, (double) nr, (double) mn, (double) mx, (double) env, (double) enu);
      fflush(stdout);
    }
    elapse += (wtm1 - wtm0);
    wtm0 = wtm1;
  }
  res->nsamp = qn;
  free(sr);

  return 0;
}
