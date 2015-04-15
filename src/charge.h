#ifndef __ETRANS_CHARGE_H
#define __ETRANS_CHARGE_H

#include <stdio.h>

#include "config.h"
#include "types.h"
#include "opt.h"

struct CHARGE_PARAM
{
  real_t chi;
  real_t lambda;
  real_t bet[4];
  real_t lam[16];
};
typedef struct CHARGE_PARAM charge_parm_t;

struct CHARGE_EQUATION
{
  int n, n0;
  real_t chi, lambda;
  real_t *d, *s, *s0; /* 3n */
  real_t *work; /* 8n */
};
typedef struct CHARGE_EQUATION charge_eq_t;

size_t charge_nbytes(int n);
charge_eq_t *mk_charge(int n, int n0, charge_parm_t *def);
int charge_init(charge_eq_t *eq, etrans_opt_t *o);
void charge_del(charge_eq_t *eq);
int charge_write(const charge_eq_t *eq, FILE *f);
int charge_read(charge_eq_t *eq, FILE *f);

void expmv(int n, real_t *bi, real_t* cr, real_t *ci, real_t *sr, real_t *si, real_t *wrk, real_t *er, real_t *ei, int *scl);


#endif /* __ETRAMS_CHARGE_H */
