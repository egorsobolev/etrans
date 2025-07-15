#ifndef __ETRANS_PBD_H
#define __ETRANS_PBD_H

#include "config.h"
#include "types.h"
#include "chain.h"
#include "charge.h"

extern charge_parm_t pbd_charge_defs;

struct PBDChain
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
  real_t omegaM2, sigma;
  real_t omegaB2, epsilon, rho;
  real_t theta0, theta;
  real_t Gamma, chi;
  real_t GammaH;
  int nheat, nmem;
  int n1;
  real_t *f0;

  real_t *runge_temp;
};
typedef struct PBDChain pbd_chain_t;

pbd_chain_t *mk_pbd(int n);

#endif /* __ETRANS_PBD_H */
