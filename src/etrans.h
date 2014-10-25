#ifndef __ETRANS_H
#define __ETRANS_H

#include <stdio.h>

#include "config.h"
#include "types.h"
#include "chain.h"
#include "charge.h"

struct EQDATA
{
	int n;
	int half;
        int h_revstep;
        real_t q_h;
        int q_nstep;
        int q_outstep;
        int x0rnd;
  /*
	real_t sv;
	real_t *d;
	real_t *s;
*/
	real_t *x0;
	/* oscilator parameters */
	real_t *x;
	real_t osc_h;
	int osc_nstep;
	int osc_outstep;
  charge_eq_t *ch;
	chain_eq_t *chain;
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
void diff_b(charge_eq_t *ch, real_t *b, real_t *u, real_t *db, real_t *d2b);

int solution_read(FILE *fn, solution_t *rs);
int solution_write(FILE *fn, solution_t *rs);
int solution_allocate(int m, int nstep, int s, solution_t *ri, solution_t *rs, double **rbuf);

#endif //__ETRANS_H
