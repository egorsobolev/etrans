#ifndef __ETRANS_PBD_H
#define __ETRANS_PBD_H

#include "config.h"
#include "types.h"
#include "chain.h"
#include "rnd.h"
#include "charge.h"

charge_parm_t pbd_charge_defs;

struct PBDChain
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
  int n, n1;
  real_t omegaM2, sigma;
  real_t omegaB2, epsilon, rho;
  real_t theta0, theta, sigmaF;
  real_t Gamma, chi;
  real_t *f0;

  real_t *runge_temp;
};
typedef struct PBDChain pbd_chain_t;

pbd_chain_t *mk_pbd(int n);

#endif /* __ETRANS_PBD_H */
