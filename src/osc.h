#ifndef __ETRANS_MORSE_H
#define __ETRANS_MORSE_H

#include "config.h"
#include "types.h"
#include "chain.h"
#include "rnd.h"

struct OSCILATOR
{
  /* inherited fields */
  autonomic_chain_step_f *step_autonomic;
  coupled_chain_step_f *step_coupled;
  autonomic_chain_eq_f *eq_autonomic;
  coupled_chain_eq_f *eq_coupled;
  chain_equilibrate_f *equilibrate;
  chain_x0_f *x0;
  chain_del_f *del;

  real_t period;
  /* own fields */
  real_t kt;
  real_t D, K, L, chi;
  real_t su, sv;
  int n, n1;
  real_t *xi;
  real_t *lambda;
  real_t *tren;
  real_t *upr;
  real_t *sig;
  real_t *runge_temp;
  rnd_stream_p rstr;
  int nsite;
};
typedef struct OSCILATOR oscilator_t;

oscilator_t *osc_init(int n, real_t temp, real_t tren, real_t upr, real_t xi, real_t lambda);

#endif /* __ETRANS_OSC_H */
