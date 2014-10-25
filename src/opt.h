#ifndef __ETRANS_OPT_H
#define __ETRANS_OPT_H

#include <argtable2.h>
#include <stdio.h>

struct ETRANS_OPTIONS
{
  struct arg_file *seqfn, *logfn, *initfn, *prmfn, *outfn;
  struct arg_lit *rst, *nxt, *help, *ver;
  struct arg_int *nsamp, *na, *nq, *ns, *nh;
  struct arg_dbl *tmax, *h, *temp, *gamma, *omega0, *chi, *mu, *xi, *hh, *drp;
  struct arg_end *end;
};
typedef struct ETRANS_OPTIONS etrans_opt_t;

void **argtable_mk(etrans_opt_t *o);
void argtable_del(void **t);
int options_write(FILE *f, int argc, char **argv);
void options_print(int argc, char **argv);

#endif /* __ETRANS_OPT_H */
