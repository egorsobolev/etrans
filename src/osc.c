#include "config.h"

#include <time.h>
#include <string.h>
#include <math.h>

#ifdef USE_MKL
# include <mkl.h>
#else
# include <gsl/gsl_cblas.h>
# include <gsl/gsl_rng.h>
# include <gsl/gsl_randist.h>
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

void f(const oscilator_t *o, const real_t *x, real_t *dx)
{
	real_t *dv;
	const real_t *v;
	int i;
	dv = dx + o->n;
	v = x + o->n;
	cblas_copy(o->n, v, 1, dx, 1);
	dv[0] = -o->upr[0] * x[0] - o->tren[0] * v[0] + o->xi[0] * x[1];
	for (i = 1; i < o->n1; ++i)
		dv[i] = -o->upr[i] * x[i] - o->tren[i] * v[i] + o->xi[i] * (x[i + 1] + x[i - 1]);
	dv[o->n1] = -o->upr[o->n1] * x[o->n1] - o->tren[o->n1] * v[o->n1] + o->xi[o->n1] * x[o->n1 - 1];
}

void f2(const oscilator_t *o, const real_t *x, real_t *dx, real_t *b2)
{
	real_t *dv;
	const real_t *v;
	int i;
	dv = dx + o->n;
	v = x + o->n;
	cblas_copy(o->n, v, 1, dx, 1);
	dv[0] = -o->upr[0] * x[0] - o->tren[0] * v[0] + o->xi[0] * x[1] - o->lambda[0] * b2[0];
	for (i = 1; i < o->n1; ++i)
		dv[i] = -o->upr[i] * x[i] - o->tren[i] * v[i] + o->xi[i] * (x[i + 1] + x[i - 1]) - o->lambda[i] * b2[i];
	dv[o->n1] = -o->upr[o->n1] * x[o->n1] - o->tren[o->n1] * v[o->n1] + o->xi[o->n1] * x[o->n1 - 1] - o->lambda[o->n1] * b2[o->n1];
}

void runge(oscilator_t *o, real_t h, real_t *x)
{
	real_t sqrth, h_2, *y, *k1, *k2;
	int i, m;
	sqrth = _sqrt(h);
	h_2 = 0.5 * h;
	m = 2 * o->n;
	y = o->runge_temp;
	k1 = y + m;
	k2 = k1 + m;

	cblas_copy(m, x, 1, y, 1);

#ifdef USE_MKL
	vRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, o->rstr, o->n, x + o->n, 0.0, 1.0);
#else
	for (i = o->n; i < m; ++i)
            x[i] = (real_t) gsl_ran_gaussian_ziggurat(o->rstr, 1.0);
#endif
	/*memset(x + o->n, 0, o->n * sizeof(float));*/

	for (i = o->n; i < m; ++i)
		x[i] = sqrth * o->sig[i - o->n] * x[i] + y[i];
	f(o, x, k1);
	cblas_axpy(m, h, k1, 1, y, 1);
	f(o, y, k2);
	cblas_axpy(m, h_2, k1, 1, x, 1);
	cblas_axpy(m, h_2, k2, 1, x, 1);
}

void runge2(oscilator_t *o, real_t h, real_t *x, real_t *b0, real_t *b1)
{
	real_t sqrth, h_2, *y, *k1, *k2;
	int i, m;
	sqrth = _sqrt(h);
	h_2 = 0.5 * h;
	m = 2 * o->n;
	y = o->runge_temp;
	k1 = y + m;
	k2 = k1 + m;

	cblas_copy(m, x, 1, y, 1);

#ifdef USE_MKL
	vRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, o->rstr, o->n, x + o->n, 0.0, 1.0);
#else
	for (i = o->n; i < m; ++i)
            x[i] = (real_t) gsl_ran_gaussian_ziggurat(o->rstr, 1.0);
#endif
	/*memset(x + o->n, 0, o->n * sizeof(float));*/

	for (i = o->n; i < m; ++i)
		x[i] = sqrth * o->sig[i - o->n] * x[i] + y[i];
	f2(o, x, k1, b0);
	cblas_axpy(m, h, k1, 1, y, 1);
	f2(o, y, k2, b1);
	cblas_axpy(m, h_2, k1, 1, x, 1);
	cblas_axpy(m, h_2, k2, 1, x, 1);
}

void osc_equilibrate(oscilator_t *osc, real_t h, real_t *x, int nstep)
{
	int i;
	for (i = 0; i < nstep; ++i)
		runge(osc, h, x);
}

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

void osc_x0_rand(oscilator_t *o, real_t *x)
{
    int i, n2;
    double sigma;
    n2 = 2 * o->n;
    sigma = _sqrt(0.5 * o->kt);
    if (o->kt == 0.0) {
        memset(x, 0, n2 * sizeof(real_t));
        return;
    }
#ifdef USE_MKL
    vRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, o->rstr, n2, x, 0.0, sigma);
#else
    for (i = 0; i < n2; ++i)
        x[i] = (real_t) gsl_ran_gaussian_ziggurat(o->rstr, sigma);
#endif
    for (i = 0; i < o->n; ++i) {
       if (o->upr[i] == 0.0)
	   x[i] = x[i + o->n] = 0.0;
       else
	   x[i] /= _sqrt(o->upr[i]);
    }
}

int osc_init(oscilator_t *o, int n, real_t temp, real_t tren, real_t upr, real_t xi, real_t lambda)
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

	o->kt = e0 / t0 * temp;
	o->n = n;
	o->n1 = n - 1;
	o->xi = (real_t *) malloc(11 * n * sizeof(real_t));
	if (!o->xi)
		return 1;
   
        o->period = 2.0 * M_PI / _sqrt(upr);

	pid = getpid();
	gettimeofday(&tm, NULL);
	seed = (int) tm.tv_usec * (int) pid;
#ifdef MPI
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	seed = seed * np + rank;
#endif
#ifdef USE_MKL
	if (vslNewStream(&o->rstr, VSL_BRNG_MT19937, seed) != VSL_STATUS_OK) {
		/*if (vslNewStream(&o->rstr, VSL_BRNG_MCG59, seed) != VSL_STATUS_OK) {*/
		free(o->xi);
		return 2;
	}
#else
        o->rstr = gsl_rng_alloc(gsl_rng_mt19937);
        if (!o->rstr) {
	   free(o->xi);
	   return 2;
	}
        gsl_rng_set(o->rstr, seed);
#endif
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
	return 0;
}

void osc_free(oscilator_t *o)
{
#ifdef USE_MKL
	vslDeleteStream(&o->rstr);
#else
        gsl_rng_free(o->rstr);
#endif
	free(o->xi);
}
