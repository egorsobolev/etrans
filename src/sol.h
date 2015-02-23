#ifndef __ETRANS_SOL_H
#define __ETRANS_SOL_H

#include <stdio.h>

#include "types.h"
#include "etrans.h"

typedef void int_f(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, double *r);
typedef long int_size_f(int n, int nstep);
typedef int int_vlen_f(int n);

struct FUCTIONAL
{
  char *name;
  int_f *f;
  int_size_f *numel;
  int_vlen_f *vlen;
};
typedef struct FUCTIONAL int_func_t;

#define MAX_FUNC (12)
int_func_t int_func[MAX_FUNC];

#define ET_MEAN (0x1)
#define ET_STD (0x2)

struct SOL_META
{
  int id;
  unsigned char flag;
  long mean_step;
  long std_step;
};

typedef struct SOL_META sol_meta_t;

struct SOL
{
  int nsamp;
  int n;
  long nstep;
  int nfunc;
  long numel;
  int nskip;
  sol_meta_t *m;
  long *smean;
  long *sstd;
  /*  int *nmean;
      int *nstd;*/
  long *imean;
  long *istd;
  int *nb;
  double *b;
  double *r;
  double **mean;
  double **std;
};
typedef struct SOL sol_t;

int sol_read_meta(FILE *f, sol_meta_t *m);
int sol_write_meta(FILE *f, int nfunc, const sol_meta_t *m);
int sol_scan_meta(FILE *f, sol_meta_t *m);
int sol_default_meta(sol_meta_t *m);
void sol_print_meta(FILE *f, int n, const sol_meta_t *m);
sol_t *mk_sol(int nfunc, sol_meta_t *m, int n, long nstep);
void sol_reset(sol_t *s);
void sol_setzero(sol_t *s);
void sol_update(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, sol_t *s);
void sol_normalize(sol_t *s, int nsamp);

#endif //__ETRANS_SOL_H
