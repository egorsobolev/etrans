#ifndef __ETRANS_INITIAL_H
#define __ETRANS_INITIAL_H

#include "types.h"

#include <stdio.h>

#ifdef MPI
#include <mpi.h>
#endif

typedef struct INITIAL_STATE initial_t;

typedef int initial_set_f(initial_t *s, int n, int n0, real_t *x);
typedef void initial_del_f(initial_t *s);
typedef size_t initial_nbytes_f(initial_t *s);

struct INITIAL_STATE
{
  initial_set_f *set;
  initial_del_f *del;
  initial_nbytes_f *nbytes;
};

struct INITIAL_1SITE
{
  initial_set_f *set;
  initial_del_f *del;
  initial_nbytes_f *nbytes;

  int n0;
};
typedef struct INITIAL_1SITE initial_1site_t;

initial_1site_t *mk_initial_1site(int n0);

struct INITIAL_1STATE
{
  initial_set_f *set;
  initial_del_f *del;
  initial_nbytes_f *nbytes;

  int n;
  real_t *x0;
};
typedef struct INITIAL_1STATE initial_1state_t;

initial_1state_t *mk_initial_1state(FILE *fn, int n, int n0);

struct INITIAL_SET
{
  initial_set_f *set;
  initial_del_f *del;
  initial_nbytes_f *nbytes;

  SFILE *f;
  int os, sl, l, n;
  real_t *x0;
};
typedef struct INITIAL_SET initial_set_t;

initial_set_t *mk_initial_set(SFILE *f, int from, int n, int n0);


#endif //__ETRANS_INITIAL_H
