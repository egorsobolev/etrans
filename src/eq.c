#include "config.h"
#include "etrans.h"
#include "osc.h"
#include "sol.h"
#include "fpset.h"

#include <rng.h>

#include <memory.h>
#include <malloc.h>
#include <stdio.h>

#include <math.h>
#include <float.h>

#ifdef MPI
#include <mpi/mpi.h>
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

size_t eq_nbytes(int n)
{
  return 23 * n * sizeof(real_t);
}

char timeunit[] = {'s', 'm', 'h'};


int eq(eqdata_t *d, real_t h, real_t tmax, int nsamp, int mg, sol_t *res)
{
  int nstep, n, n1, nn, i, j, ln, k, q, qn, l;
  real_t a, p2, normp, nr, sqrth, normA, t;
  real_t *u, *v, *b2l, *b2r, *b2t;
  real_t *dsr, *d2sr, *x3, *x2, *x1, *xt, *u3, *u2, *u1, *ut, *g1, *g2;
  real_t *di, *cr, *ci;
  int rank, err, progress, tu;
  double elapse, wtm0, wtm1, dwtm;

#ifdef MPI
  MPI_Status status;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  rank = 0;
#endif


  qn = nsamp;
  n = d->n;
  n1 = n - 1;
  nn = 2*n;
  nstep = d->nstep;

  sqrth = _sqrt(h) * d->chain->sigF;

  di = (real_t *) malloc(eq_nbytes(n));
  if (!di) return -1;

  cr = di + n;
  ci = cr + n;
  x3 = ci + n;
  u3 = x3 + nn;
  x2 = u3 + nn;
  u2 = x2 + nn;
  x1 = u2 + nn;
  u1 = x1 + nn;
  g1 = u1 + nn;
  dsr = g1 + nn;
  d2sr = dsr + nn;
  b2l = d2sr + nn;
  b2r = b2l + n;

  sol_setzero(res);
   
  elapse = 0;
  wtm0 = walltime();
  for (q = 0; q < qn; q++) {

    d->chain->x0(d->chain, u3);

    if (d->s->set(d->s, n, d->ch->n0, x3))
      break;

    if (!rank) {
      printf("%8d ", q + 1);
      fflush(stdout);
    }

    p2 = 0.0;
    for (k = 0; k < n; ++k) {
      a = x3[k] * x3[k] + x3[k+n] * x3[k+n];
      b2l[k] = a;
      p2 += a;
    }
    nr = _fabs(p2 - 1.0);

    diff_b(d->ch, x3, dsr, d2sr);
    sol_reset(res);
    sol_update(d, x3, b2l, dsr, d2sr, res);

    j = d->nskip - 1;
    ln = 0;
    l = 0;
    normA = 0.0;
    t = d->ntm;

    /* cycle by classical steps */
    i = 0;
    while (i < nstep) {
      xt = x1; x1 = x2; x2 = x3; x3 = xt;
      ut = u1; u1 = u2; u2 = u3; u3 = ut;

      /* 2o2s1g algorithm */
      /* HS Greenside and E Helfand, The Bell System Technical Journal 60 (1981) 1927-1940 */
      rng_gaussian(n, u3 + n, 1.0);
      for (k = 0; k < n; ++k) {
	u3[k] = u2[k];
	u3[k+n] = sqrth * u3[k+n] + u2[k+n];
	g1[k] = u3[k+n]; /* du = v */
      }
      d->chain->dv_hst(d->chain, u3, g1+n, b2l);
      for (k = 0; k < nn; ++k) {
	g1[k] *= h;
	u3[k] += 0.5 * g1[k];
	g1[k] += u2[k];
      }
      for (k = 0; k < n; ++k)
	u3[k] += 0.5 * h * g1[k+n];

      if (mg == 2 && i) {
	/* Magnus 4rd order */
	/* S Blanes et al./ Physics Reports 470 (2009) 151-238 */
	/* h * (A1 + 4A2 + A3) / 6 - h^2 [A1, A3] / 12 */
	for (k = 0; k < n; ++k) {
	  di[k] = 2 * h * (d->ch->d[k] + d->ch->chi * (u1[k] + 4 * u2[k] + u3[k]) / 6);
	}
	for (k = 0; k < n1; ++k) {
	  ci[k] = 2 * h * d->ch->s[k];
	  cr[k] = h * h * d->ch->s[k] * d->ch->chi * ((u1[k] - u3[k]) - (u1[k+1] - u3[k+1])) / 3;
	}
      } else if (mg == 1 || (mg == 2 && ~i)) {
	/* Magnus 2nd order */
	/* S Blanes et al./ Physics Reports 470 (2009) 151-238 */
	/* h * (A1 + A2) / 2 */
	for (k = 0; k < n; ++k) {
	  di[k] = h * (d->ch->d[k] + 0.5 * d->ch->chi * (u2[k] + u3[k]));
	}
	for (k = 0; k < n1; ++k) {
	  ci[k] = h * d->ch->s[k];
	  cr[k] = 0.0;
	}
	for (k = 0; k < nn; ++k)
	  x1[k] = x2[k];
      } else {
        /* 1st order sheme */
        /* guess next u by current point */
	for (k = 0; k < n; ++k) {
	  di[k] = h * (d->ch->d[k] + d->ch->chi * (u2[k] + 0.5 * h * u2[k+n]));
	}
	for (k = 0; k < n1; ++k) {
	  ci[k] = h * d->ch->s[k];
	  cr[k] = 0.0;
	}
	for (k = 0; k < nn; ++k)
	  x1[k] = x2[k];
      }
      /* exponential matrix times vector
	 G Gallopoulos and Y Saad, SIAM Journal on Scientific and Statistical Computing 13 (1992) 1236-1264

	 NB: Here dsr(2n) and d2sr(2n) are used as workspace and destroyed */
      expmv(n, di, cr, ci, x1, x1+n, dsr, x3, x3+n);

      if (*dsr > normA)
	normA = *dsr;

      /* calculate I1 = sum(b(k)^2) and move b(k) back to the trajectory */
      p2 = 0.0;
      for (k = 0; k < n; ++k) {
	a = x3[k] * x3[k] + x3[n+k] * x3[n+k];
	b2r[k] = a;
	p2 += a;
      }
      normp = _fabs(p2 - 1.0);
      if (normp > nr)
	nr = normp;
      a = _sqrt(p2);

      if (normp > d->nthr || t <= 0.0) {
	for (k = 0; k < n; ++k) {
	  x3[k] /= a;
	  x3[k+n] /= a;
	  b2r[k] /= p2;
	}
	if (t <= 0.0)
	  t = d->ntm;
      }
      /* 2o2s1g algorithm continue 
	 HS Greenside and E Helfand, The Bell System Technical Journal 60 (1981) 1927-1940 

         NB: Here u1(n:2n) are used as workspace and destroyed */
      d->chain->dv_hst(d->chain, g1, u1+n, b2r);
      for (k = n; k < nn; ++k)
	u3[k] += 0.5 * h * u1[k];

      b2t = b2l;
      b2l = b2r;
      b2r = b2t;

      if (!j) {
	diff_b(d->ch, x3, dsr, d2sr);
	sol_update(d, x3, b2l, dsr, d2sr, res);
	j = d->nskip;
	++i;
      }

      if (!rank) {
	progress = (int) (22.0 * (1.0 + i - j / (double) d->nskip) / (double) nstep + 0.5);
	while (l < progress) {
	  printf(".");
	  fflush(stdout);
	  l++;
	}
      }

      --j;
      --ln;
      t -= h;
    }
    if (fpset_write(4*n, x3)) {
      free(di);
      return -2;
    }
    wtm1 = walltime();
    dwtm = wtm1 - wtm0;
    tu = 0;
    while (tu < 2 && dwtm > 60) {
      dwtm /= 60;
      tu++;
    }
    /*               1         2         3         4         5         6         7         8
    /*                12345678901234567890123456789012345678901234567890123456789012345678901234567890 */
    /* printf("     # ===== progress =====      t,s max|P2-1|    min(h)     P(1)     P(n)\n"); */
    if (!rank) {
      printf(" %8.3g%c%10.2g%9.2g%9.1e%9.1e\n", dwtm, timeunit[tu], (double) nr, (double) normA, (double) b2l[0], (double) b2l[n-1]);
      fflush(stdout);
    }
    elapse += (wtm1 - wtm0);
    wtm0 = wtm1;
  }
  res->nsamp = q;
  free(di);

  return 0;
}
