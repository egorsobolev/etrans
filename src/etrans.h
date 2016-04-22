#ifndef __ETRANS_H
#define __ETRANS_H

#include <stdio.h>

#include "config.h"
#include "types.h"
#include "chain.h"
#include "charge.h"
#include "initial.h"

struct EQDATA
{
  int n;
  long nstep;
  real_t h;
  int nskip;

  real_t ntm;
  real_t nthr;

  charge_eq_t *ch;
  chain_eq_t *chain;
  charge_parm_t *charge_parm;
  initial_t *s;
};
typedef struct EQDATA eqdata_t;

int readparm(FILE *f, charge_parm_t *p);
int seqscan(FILE *f, char **seq);
int seqdna(charge_parm_t *p, int n, const char *seq, real_t *d, real_t *s);

/*int b0_read(FILE *f, int n, real_t *x);*/
void diff_b(charge_eq_t *ch, real_t *b, real_t *db, real_t *d2b);

size_t eq_nbytes(int n);

#endif //__ETRANS_H
