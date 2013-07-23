#include "config.h"
#include "etrans.h"

#include <memory.h>
#include <malloc.h>
#include <stdio.h>

#include <math.h>
#include <float.h>

#ifdef USE_MKL
# include <mkl.h>
#else
# include <gsl/gsl_cblas.h>
# include <gsl/gsl_rng.h>
#endif

#ifdef MPI
#include <mpi.h>
#define walltime() MPI_Wtime()
#else
#include <time.h>
#define walltime()	((double) clock() / CLOCKS_PER_SEC)
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

int eq(eqdata_t *d, real_t h, real_t tmax, int nsamp, solution_t *r)
{
	int nstep, nskip, osc_n1, n, i, j, jn, ln, k, l, m, q, qn, tn, t, neigs, info, ps;
	real_t p2, w, abstol, a, b, t0, sn, cs, h0, hmax, normp, mn, mx, nr, osc_period, osc_h;
	real_t *u, *v, *s0, *d0, *ds, *dk, *s0r, *s0i, *tr, *ti, *work, *wa, *b2l, *b2r, *b2t, *x;
        real_t *sr, *si, *dsr, *dsi, *d2sr, *d2si;
	int *iwork, rank;
	double elapse, wtm0, wtm1, c, dc, d2c;

#ifdef MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = 0;
#endif

	abstol = 2.0 * clapack_lamch('S');
	osc_n1 = d->osc_nstep - 1;

	qn = r->nsamp + nsamp;
	n = d->n;
        r->nstep = (int) ceil(tmax / (h * d->q_outstep));
	nstep = r->nstep * d->q_outstep;
	r->n = n;

	d0 = (real_t *) malloc(24 * n * sizeof(real_t) + 5 * n * sizeof(int));
	if (!d0) return -1;

	s0 = d0 + n;
	sr = s0 + n;
	si = sr + n;
	s0r = si + n;
	s0i = s0r + n;
        dsr = s0i + n;
        dsi = dsr + n;
        d2sr = dsi + n;
        d2si = d2sr + n;
	dk = d2si + n;
	ds = dk + n;
	tr = ds + n;
	ti = tr + n;
	wa = ti + n;
	x = wa + n;
	b2l = x + 2 * n;
	b2r = b2l + n;

	work = b2r + n;
	iwork = (int *) (work + 5 * n);

	u = x;
	v = u + n;
   
        osc_h = d->osc.period / 1000.0;
   
	elapse = 0;
	wtm0 = walltime();
	for (q = r->nsamp; q < qn; q++) {

	        if (d->x0rnd) {
	            osc_x0_rand(&d->osc, d->x);
                } else {
#ifdef USE_MKL
		    viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, d->osc.rstr, 1, &ps, 0, 1000);
#else
		    ps = gsl_rng_uniform_int(d->osc.rstr, 1000);
#endif
		    osc_equilibrate(&d->osc, osc_h, d->x, ps);
                }
		cblas_copy(2 * n, d->x, 1, x, 1);
		cblas_axpy(2 * n, 1.0, d->x0 + 2 * n, 1, x, 1); 
		cblas_copy(2 * n, d->x0, 1, sr, 1);
		p2 = 1.0;

		diff_b(d, sr, x, dsr, d2sr);
		c = 0.0;
		dc = 0.0;
		d2c = 0.0;
		for (k = 0; k < n; ++k) {
			a = sr[k] * sr[k] + si[k] * si[k];
			b2l[k] = a;
			b = (real_t) k - d->half;
		        b = b * b;
			c += (double) a * b;
		        b *= 2.0;
			dc += (double) (sr[k] * dsr[k] + si[k] * dsi[k]) * b;
			d2c += (double) (dsr[k] * dsr[k] + dsi[k] * dsi[k] + sr[k] * d2sr[k] + si[k] * d2si[k]) * b;
			
		        r->p[k] += (double) a;
			r->p2[k] += (double) a * a;

			r->u[k] += (double) u[k];
			r->u2[k] += (double) u[k] * u[k];

			p2 += a;
		}
		r->c[0] += c;
		r->dc[0] += dc;
		r->d2c[0] += d2c;
		r->c2[0] += c * c;
		r->dc2[0] += dc * dc;
		r->d2c2[0] += d2c * d2c;

		nr = 0.0;
		mn = REAL_MAX;
		mx = 0.0;

		jn = d->osc_outstep;
	        tn = d->q_outstep;
	        ln = 0;
		j = 0;
	        t = 0;
		/**********************/

		/* cycle by classical steps */
		for (i = 0; i < nstep; ++i) {

		        if (!ln) {
				/* d0 = d + u(j) */
				cblas_copy(n, d->d, 1, d0, 1);
				cblas_axpy(n, d->sv, u, 1, d0, 1);
				cblas_copy(n, d0, 1, ds, 1);
				cblas_axpy(n, d->h_revstep * h * d->sv, v, 1, ds, 1);

				clapack_stevx('N', 'I', n, d0, d->s, 0.0, 1.0, 1, 1, abstol, &neigs, wa, NULL, n, work, iwork, NULL, &info);
				a = wa[0];
				clapack_stevx('N', 'I', n, ds, d->s, 0.0, 1.0, 1, 1, abstol, &neigs, wa, NULL, n, work, iwork, NULL, &info);
				b = wa[0];

				w = a < b ? a : b;
				for (k = 0; k < n; ++k) {
					d0[k] -= w;
					ds[k] -= w;
				}
				clapack_stevx('N', 'I', n, d0, d->s, 0.0, 1.0, n, n, abstol, &neigs, wa, NULL, n, work, iwork, NULL, &info);
				a = wa[0];
				clapack_stevx('N', 'I', n, ds, d->s, 0.0, 1.0, n, n, abstol, &neigs, wa, NULL, n, work, iwork, NULL, &info);
				b = wa[0];

				hmax = 0.75 * M_PI / (a > b ? a : b);
				if (hmax > d->osc_h) hmax = d->osc_h;
				nskip = h > hmax ? (int) ceil(h / hmax) : 1;
				h0 = h / nskip;
				if (h0 > mx) mx = h0;
				if (h0 < mn) mn = h0;

				a = w * h0;
				sn = 0.5 * _sin(a);
				cs = 0.5 * _cos(a);

				cblas_copy(n - 1, d->s, 1, s0, 1);
				cblas_scal(n - 1, h0, s0, 1);
			   
			        ln = d->h_revstep;
			}
			for (m = 0; m < nskip; ++m) {
				/* ds = h * (d0 + dt * v(j)) */
				cblas_copy(n, d->d, 1, ds, 1);
				cblas_axpy(n, d->sv, u, 1, ds, 1);
				/*cblas_scopy(n, d0, 1, ds, 1);*/
				cblas_axpy(n, (m + 0.5) * h0 * d->sv, v, 1, ds, 1);
				cblas_scal(n, h0, ds, 1);

				cblas_copy(n, sr, 1, s0r, 1);
				cblas_copy(n, si, 1, s0i, 1);

				cblas_scal(n, a0, s0r, 1);
				cblas_scal(n, a0, s0i, 1);
				for (l = 0; l < 14; ++l) {
					for (k = 0; k < n; ++k)
						dk[k] = ds[k] - thi[l];
					thomas1(n, dk, -thr[l], s0, si, sr, ti, tr, work);
					/* s0r = s0r + ar[k] * tr - ai[k] * ti */
					cblas_axpy(n,  ar[l], tr, 1, s0r, 1);
					cblas_axpy(n, -ai[l], ti, 1, s0r, 1);
					/* s0i = s0i + ai[k] * tr + ar[k] * ti */
					cblas_axpy(n,  ai[l], tr, 1, s0i, 1);
					cblas_axpy(n,  ar[l], ti, 1, s0i, 1);
				}
				cblas_copy(n, s0r, 1, sr, 1);
				cblas_copy(n, s0i, 1, si, 1);
				/* sr = cs * s0r + sn * s0i */
				cblas_scal(n, cs, sr, 1);
				cblas_axpy(n, sn, s0i, 1, sr, 1);
				/* si = cs * s0i - sn * s0r */
				cblas_scal(n,  cs, si, 1);
				cblas_axpy(n, -sn, s0r, 1, si, 1);

				p2 = cblas_dot(2 * n, sr, 1, sr, 1);
				normp = _fabs(p2 - 1.0);
				if (normp > nr)
					nr = normp;

				p2 = _sqrt(p2);
				a = 1.0 / p2;
				cblas_scal(2 * n, a, sr, 1);

/*				t0 += h0;*/
			}
			for (k = 0; k < n; ++k) {
				b2r[k] = sr[k] * sr[k] + si[k] * si[k];
			}
			runge2(&d->osc, h, x, b2l, b2r);
			b2t = b2l;
			b2l = b2r;
			b2r = b2t;
		     
		        --tn;
		        if (!tn) {
                           ++t;
			   
			   diff_b(d, sr, x, dsr, d2sr);
			   c = 0.0;
			   dc = 0.0;
			   d2c = 0.0;
			   for (k = 0; k < n; ++k) {
			      a = b2l[k];
			      b = (real_t) k - d->half;
			      b = b * b;
			      c += (double) a * b;
			      b *= 2.0;
			      dc += (double) (sr[k] * dsr[k] + si[k] * dsi[k]) * b;
			      d2c += (double) (dsr[k] * dsr[k] + dsi[k] * dsi[k] + sr[k] * d2sr[k] + si[k] * d2si[k]) * b;
			      r->p[k + t * n] += (double) a;
			      r->p2[k + t * n] += (double) a * a;
			   } 
			   r->c[t] += c;
			   r->dc[t] += dc;
			   r->d2c[t] += d2c;
			   r->c2[t] += c * c;
			   r->dc2[t] += dc * dc;
			   r->d2c2[t] += d2c * d2c;
			   
			   tn = d->q_outstep;
			}
		   
			--jn;
			if (!jn) {
				++j;
		                for (k = 0; k < n; ++k) {
					r->u[j * n + k] += (double) u[k];
					r->u2[j * n + k] += (double) u[k] * u[k];
				}
				jn = d->osc_outstep;
			}		   
		        --ln;
		}
		wtm1 = walltime();
		if (!rank) {
		        c = _sqrt(0.5 * n * d->osc.kt);
		        a = 2.0 * M_PI * cblas_nrm2(n, d->x, 1) / (c *  d->osc.period);
		        b = cblas_nrm2(n, d->x + n, 1) / c;
			printf("%6d %8.2f %9.2g %9.6g %9.6g %9.6g %9.6g\n", q + 1, wtm1 - wtm0, (double) nr, (double) mn, (double) mx, (double) a, (double) b);
			fflush(stdout);
		}
		elapse += (wtm1 - wtm0);
		wtm0 = wtm1;
	}
	r->nsamp = qn;
	free(d0);
        /*vslDeleteStream(&pstr);*/

	return 0;
}
