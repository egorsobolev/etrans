#ifndef __ETRANS_OSC_H
#define __ETRANS_OSC_H

#include "config.h"
#include "types.h"
#include "chain.h"
#include "rnd.h"
#include "charge.h"

charge_parm_t osc_charge_defs;

struct OSCILATOR
{
  /* inherited fields */
  chain_dv_f *dv;
  chain_dv_hst_f *dv_hst;
  chain_equilibrate_f *equilibrate;
  en_potential_f *en_potential;
  chain_x0_f *x0;
  chain_del_f *del;
  chain_init_f *init;
  chain_write_f *write;
  chain_read_f *read;
  chain_nbytes_f *nbytes;

  int n;
  real_t h;
  real_t sigF;

  /* own fields */
  real_t kt;
  real_t D, K, L, mu;
  real_t su, sv;
  int n1;
};
typedef struct OSCILATOR oscilator_t;

oscilator_t *mk_osc(int n);

#endif /* __ETRANS_OSC_H */
