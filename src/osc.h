#ifndef __ETRANS_MORSE_H
#define __ETRANS_MORSE_H

#include "config.h"
#include "types.h"
#include "chain.h"
#include "rnd.h"
#include "charge.h"

charge_parm_t osc_charge_defs;

struct OSCILATOR
{
  /* inherited fields */
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

  /* own fields */
  real_t kt;
  real_t D, K, L, mu, sigF;
  real_t su, sv;
  int n, n1;

  real_t *runge_temp;
};
typedef struct OSCILATOR oscilator_t;

oscilator_t *mk_osc(int n);

#endif /* __ETRANS_OSC_H */
