#ifndef __ETRANS_H
#define __ETRANS_H

#include <stdio.h>

#include "config.h"

#ifdef USE_MKL
# include <mkl.h>
typedef VSLStreamStatePtr rng_stream;
#else
# include <gsl/gsl_rng.h>
typedef gsl_rng *rng_stream;
#endif

#ifdef USE_SINGLE
typedef float real_t;

# define _sqrt(X) sqrtf(X)
# define _sin(X) sinf(X)
# define _cos(X) cosf(X)
# define _fabs(X) fabsf(X)
# define cblas_copy cblas_scopy
# define cblas_scal cblas_sscal
# define cblas_dot cblas_sdot
# define cblas_nrm2 cblas_snrm2
# define cblas_axpy cblas_saxpy
# define clapack_lamch clapack_slamch
# define clapack_stevx clapack_sstevx
# define vRngGaussian vsRngGaussian
# define thomas1 thomas1s
#else
typedef double real_t;

# define _sqrt(X) sqrt(X)
# define _sin(X) sin(X)
# define _cos(X) cos(X)
# define _fabs(X) fabs(X)
# define cblas_copy cblas_dcopy
# define cblas_scal cblas_dscal
# define cblas_dot cblas_ddot
# define cblas_nrm2 cblas_dnrm2
# define cblas_axpy cblas_daxpy
# define clapack_lamch clapack_dlamch
# define clapack_stevx clapack_dstevx
# define vRngGaussian vdRngGaussian
# define thomas1 thomas1d
#endif

struct OSCILATOR
{
	real_t kt, period;
	int n, n1;
	real_t *xi;
	real_t *lambda;
	real_t *tren;
	real_t *upr;
	real_t *sig;
	real_t *runge_temp;
	rng_stream rstr;
	int nsite;
};
typedef struct OSCILATOR oscilator_t;

struct EQDATA
{
	int n;
	int half;
        int h_revstep;
        real_t q_h;
        int q_nstep;
        int q_outstep;
        int x0rnd;
	real_t sv;
	real_t *d;
	real_t *s;
	real_t *x0;
	/* oscilator parameters */
	real_t *x;
	real_t osc_h;
	int osc_nstep;
	int osc_outstep;
	oscilator_t osc;
};
typedef struct EQDATA eqdata_t;

struct SOLUTION
{
	int n;
	int nstep;
	int nsamp;
	float h;
        int pos;
        int uos;
	double *p;
	double *c;
	double *dc;
	double *d2c;
	double *u;
	double *p2;
	double *c2;
	double *dc2;
	double *d2c2;
	double *u2;
/*        float hw;
        int *hist;*/
};
typedef struct SOLUTION solution_t;

int readparm(FILE *f);
int seqscan(FILE *f, char **seq);
int seqdna(int n, const char *seq, real_t *d, real_t *s);

int eq(eqdata_t *d, real_t h, real_t tmax, int nsamp, solution_t *r);

int b0_read(FILE *f, int n, real_t *x);
void diff_b(eqdata_t *d, real_t *b, real_t *u, real_t *db, real_t *d2b);

int osc_init(oscilator_t *o, int n, real_t temp, real_t tren, real_t upr, real_t xi, real_t lambda);
void osc_equilibrate(oscilator_t *osc, real_t h, real_t *x, int nstep);
void osc_integrate(oscilator_t *osc, real_t h, real_t *x, int nskip, int nout, real_t *u, real_t *v);
void osc_free(oscilator_t *o);
void osc_x0_rand(oscilator_t *o, real_t *x);

void runge(oscilator_t *o, real_t h, real_t *x);
void runge2(oscilator_t *o, real_t h, real_t *x, real_t *b0, real_t *b1);
void f(const oscilator_t *o, const real_t *x, real_t *dx);
void f2(const oscilator_t *o, const real_t *x, real_t *dx, real_t *b2);

int solution_read(FILE *fn, solution_t *rs);
int solution_write(FILE *fn, solution_t *rs);
int solution_allocate(int m, int nstep, int s, solution_t *ri, solution_t *rs, double **rbuf);

#endif //__ETRANS_H
