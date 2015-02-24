#ifndef __ETRANS_CHAIN_H
#define __ETRANS_CHAIN_H

#include "config.h"
#include "types.h"

#include "opt.h"

typedef struct CHAIN_EQUATION chain_eq_t;

typedef void autonomic_chain_step_f(chain_eq_t *o, real_t h, real_t *x);
typedef void coupled_chain_step_f(chain_eq_t *o, real_t h, real_t *x, real_t *b0, real_t *b1);
typedef void autonomic_chain_eq_f(const chain_eq_t *o, const real_t *x, real_t *dx);
typedef void coupled_chain_eq_f(const chain_eq_t *o, const real_t *x, real_t *dx, real_t *b2);

typedef void chain_equilibrate_f(chain_eq_t *chain);
typedef void chain_del_f(chain_eq_t *o);
typedef void chain_x0_f(chain_eq_t *o, real_t *x);

typedef int chain_init_f(chain_eq_t *c, etrans_opt_t *o);

typedef real_t en_potential_f(const chain_eq_t *o, const real_t *u);

typedef int chain_write_f(const chain_eq_t *c, FILE *f);
typedef int chain_read_f(chain_eq_t *c, FILE *f);

typedef size_t chain_nbytes_f(int n);

struct CHAIN_EQUATION
{
  autonomic_chain_step_f *step_autonomic;
  coupled_chain_step_f *step_coupled;
  autonomic_chain_eq_f *eq_autonomic;
  coupled_chain_eq_f *eq_coupled;
  chain_equilibrate_f *equilibrate;
  en_potential_f *en_potential;
  chain_x0_f *x0;
  chain_del_f *del;
  chain_init_f *init;
  chain_write_f *write;
  chain_read_f *read;
  chain_nbytes_f *nbytes;

  real_t h;
};


#endif /*  __ETRANS_CHAIN_H */
