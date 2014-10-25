#ifndef __ETRANS_CHAIN_H
#define __ETRANS_CHAIN_H

#include "config.h"
#include "types.h"

typedef struct CHAIN_EQUATION chain_eq_t;

typedef void autonomic_chain_step_f(chain_eq_t *o, real_t h, real_t *x);
typedef void coupled_chain_step_f(chain_eq_t *o, real_t h, real_t *x, real_t *b0, real_t *b1);
typedef void autonomic_chain_eq_f(const chain_eq_t *o, const real_t *x, real_t *dx);
typedef void coupled_chain_eq_f(const chain_eq_t *o, const real_t *x, real_t *dx, real_t *b2);

typedef void chain_equilibrate_f(chain_eq_t *osc, real_t h, real_t *x, int nstep, int rstep);
typedef void chain_del_f(chain_eq_t *o);
typedef void chain_x0_f(chain_eq_t *o, real_t *x);

struct CHAIN_EQUATION
{
  autonomic_chain_step_f *step_autonomic;
  coupled_chain_step_f *step_coupled;
  autonomic_chain_eq_f *eq_autonomic;
  coupled_chain_eq_f *eq_coupled;
  chain_equilibrate_f *equilibrate;
  chain_x0_f *x0;
  chain_del_f *del;

  real_t period;
};


#endif /*  __ETRANS_CHAIN_H */
