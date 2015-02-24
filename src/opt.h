#ifndef __ETRANS_OPT_H
#define __ETRANS_OPT_H

#include <argtable2.h>
#include <stdio.h>

struct ETRANS_OPTIONS
{
  struct arg_file *seqfn, *logfn, *initfn, *prmfn, *outfn, *outsh, *lpfn;
  struct arg_lit *rst, *help, *ver;
  struct arg_int *nsamp, *na, *n0, *no;
  struct arg_dbl *tmax, *h, *temp, *gamma, *omegaM2, *chi, *mu, *omegaB2, *drp, *sigma, *rho, *epsilon, *lambda;
  struct arg_str *mdl, *init;
  struct arg_end *end;
};
typedef struct ETRANS_OPTIONS etrans_opt_t;

void **argtable_mk(etrans_opt_t *o);
void argtable_del(void **t);
int options_write(FILE *f, int argc, char **argv);
void options_print(int argc, char **argv);

#endif /* __ETRANS_OPT_H */
