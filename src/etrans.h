#ifndef __ETRANS_H
#define __ETRANS_H

#include <mkl.h>
#include <stdio.h>

#define REFERENCE_OPS 2e6

struct OSCILATOR
{
	float kt, period;
	int n, n1;
	float *xi;
	float *lambda;
	float *tren;
	float *upr;
	float *sig;
	float *runge_temp;
	VSLStreamStatePtr rstr;
	int nsite;
};
typedef struct OSCILATOR oscilator_t;

struct EQDATA
{
	int n;
	int half;
        int h_revstep;
        int q_outstep;
        int x0rnd;
	float sv;
	float *d;
	float *s;
	float *x0;
	/* oscilator parameters */
	float *x;
	float osc_h;
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
	int o_nstep;
	float o_h;
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
int seqdna(int n, const char *seq, float *d, float *s);

int eq(eqdata_t *d, float h, float tmax, int nsamp, solution_t *r);

int b0_read(FILE *f, int n, float *x);
void diff_b(eqdata_t *d, float *b, float *u, float *db, float *d2b);

int osc_init(oscilator_t *o, int n, float temp, float tren, float upr, float xi, float lambda);
void osc_equilibrate(oscilator_t *osc, float h, float *x, int nstep);
void osc_integrate(oscilator_t *osc, float h, float *x, int nskip, int nout, float *u, float *v);
void osc_free(oscilator_t *o);
void osc_x0_rand(oscilator_t *o, float *x);

void runge(oscilator_t *o, float h, float *x);
void runge2(oscilator_t *o, float h, float *x, float *b0, float *b1);
void f(const oscilator_t *o, const float *x, float *dx);
void f2(const oscilator_t *o, const float *x, float *dx, float *b2);

int solution_read(FILE *fn, solution_t *rs);
int solution_write(FILE *fn, solution_t *rs);
int solution_allocate(int m, int nstep, int s, solution_t *ri, solution_t *rs, double **rbuf);

#endif //__ETRANS_H
