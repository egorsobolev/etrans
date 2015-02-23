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

charge_eq_t *mk_charge(int n, int n0, charge_parm_t *def);
int charge_init(charge_eq_t *eq, etrans_opt_t *o);
void charge_del(charge_eq_t *eq);

real_t charge_expM(charge_eq_t *eq, int m, real_t h0, real_t sn, real_t cs, const real_t *u, const real_t *v, real_t *sr, real_t *si);
void charge_step_adj(charge_eq_t *eq, real_t h, int revstep, const real_t *u, const real_t *v, const real_t *sr, const real_t *si,
                        int *nskip, real_t *h0, real_t *sn, real_t *cs, real_t *s0);

int charge_write(const charge_eq_t *eq, FILE *f);
int charge_read(charge_eq_t *eq, FILE *f);

#endif /* __ETRAMS_CHARGE_H */
