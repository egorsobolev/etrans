#ifndef __ETRANS_CHARGE_H
#define __ETRANS_CHARGE_H

#include "config.h"
#include "types.h"

struct CHARGE_EQUATION
{
  int n;
  real_t chi;
  real_t *d, *s, *s0; /* 3n */
  real_t *work; /* 8n */
  int h_revstep;
};
typedef struct CHARGE_EQUATION charge_eq_t;

int charge_init(charge_eq_t *eq, int n, real_t chi);
void charge_del(charge_eq_t *eq);
real_t charge_expM(charge_eq_t *eq, int m, real_t h0, real_t sn, real_t cs, const real_t *u, const real_t *v, real_t *sr, real_t *si);
void charge_step_adj(charge_eq_t *eq, real_t h, const real_t *u, const real_t *v, const real_t *sr, const real_t *si,
                        int *nskip, real_t *h0, real_t *sn, real_t *cs, real_t *s0);

#endif /* __ETRAMS_CHARGE_H */
