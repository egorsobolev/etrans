#include "initial.h"
#include "rnd.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef USE_MKL
# include <mkl.h>
#else
# include <gsl/gsl_cblas.h>
#endif

void initial_del(initial_t *state)
{
  free(state);
}

int site_set(initial_t *state, int n, int n0, real_t *x)
{
  initial_1site_t *s = (initial_1site_t *) state;
  real_t a;
  real_t *y;

  y = x + n;

  memset(x, 0, 2 * n * sizeof(real_t));
  rnd_uniform(1, &a, 2.0 * M_PI);
  x[n0] = cos(a);
  y[n0] = sin(a);

  return 0;
}

size_t site_nbytes_impl()
{
  return sizeof(initial_1site_t);
}

size_t site_nbytes(initial_t *state)
{
  return site_nbytes_impl();
}

initial_1site_t *mk_initial_1site(int n0)
{
  initial_1site_t *s;
  s = (initial_1site_t *) malloc(site_nbytes_impl());
  if (!s)
    return NULL;
  s->n0 = n0;
  s->set = &site_set;
  s->del = &initial_del;
  s->nbytes = &site_nbytes;
  
  return s;
}


int state_set(initial_t *state, int n, int n0, real_t *x)
{
  initial_1state_t *s = (initial_1state_t *) state;
  int nn = 2 * n;

  cblas_axpy(nn, 1.0, s->x0 + nn, 1, x + nn, 1);
  cblas_copy(nn, s->x0, 1, x, 1);

  return 0;
}

size_t state_nbytes_impl(int n)
{
  return sizeof(initial_1state_t) + 4 * n * sizeof(real_t);
}

size_t state_nbytes(initial_t *state)
{
  initial_1state_t *s = (initial_1state_t *) state;
  return state_nbytes_impl(s->n);
}

initial_1state_t *mk_initial_1state(FILE *f, int n, int n0)
{
  initial_1state_t *s;
  int m, m0, sr, sl, os, l, i;
  real_t *xi;

  s = (initial_1state_t *) malloc(state_nbytes_impl(n));
  if (!s)
    return NULL;
  s->x0 = (real_t *) (s + 1);
  s->set = &state_set;
  s->del = &initial_del;
  s->nbytes = &state_nbytes;

  memset(s->x0, 0, 4 * n * sizeof(real_t));
  s->n = n;

  if (fread(&m, sizeof(int), 1, f) != 1)
    goto err;
  if (fread(&m0, sizeof(int), 1, f) != 1)
    goto err;

  os = n0 - m0;
  l = n - m;

  sr = os < l ? 0 : os - l;
  if (os > 0) {
    sl = 0;
  } else {
    sl = -os;
    os = 0;
  }
  sr += sl;
  l = m - sr;
  xi = s->x0 + os;

  if (fseek(f, sl, SEEK_CUR))
    goto err;
  for (i = 0; i < 3; i++) {
    if (fread(xi, sizeof(real_t), l, f) != l)
      goto err;
    if (fseek(f, sr, SEEK_CUR))
      goto err;
    xi += n;
  }
  if (fread(xi, sizeof(real_t), l, f) != l)
    goto err;

  return s;
 err:
  free(s);
  return NULL;
}


int set_set(initial_t *state, int n, int n0, real_t *x)
{
  initial_set_t *s = (initial_set_t *) state;
#ifdef MPI
  if (MPI_File_read_shared(*s->f, s->x0, 4 * s->n, MPI_ET_REAL, MPI_STATUS_IGNORE))
    return -1;
#else
  if (fread(s->x0, sizeof(real_t), 4 * s->n, s->f) != (4 * s->n))
    return -1;
#endif
  memset(x, 0, 2 * n * sizeof(real_t));

  cblas_copy(s->l, s->x0+s->sl, 1, x+s->os, 1);
  cblas_copy(s->l, s->x0+s->n+s->sl, 1, x+n+s->os, 1);

  cblas_copy(s->l, s->x0+2*s->n+s->sl, 1, x+2*n+s->os, 1);
  cblas_copy(s->l, s->x0+3*s->n+s->sl, 1, x+3*n+s->os, 1);

  return 0;
}

size_t set_nbytes_impl(int n)
{
  return sizeof(initial_set_t) + 4 * n * sizeof(real_t);
}

size_t set_nbytes(initial_t *state)
{
  initial_set_t *s = (initial_set_t *) state;
  return set_nbytes_impl(s->n);
}

initial_set_t *mk_initial_set(SFILE *f, int from, int n, int n0)
{
  initial_set_t *s;
  int b[2], os, sr, l;

#ifdef MPI
  if (MPI_File_read_at_all(*f, 0, b, 2, MPI_ET_REAL, MPI_STATUS_IGNORE))
    return NULL;
  if (MPI_File_seek_shared(*f, 2 * sizeof(int) + 4 * b[0] * from, MPI_SEEK_SET))
    return NULL;
#else
  if (fread(&b, sizeof(int), 2, f) != 2)
    return NULL;
  if (fseek(f, 4 * b[0] * from, SEEK_CUR))
    return NULL;
#endif

  s = (initial_set_t *) malloc(set_nbytes_impl(b[0]));
  if (!s)
    return NULL;
  s->x0 = (real_t *) (s + 1);

  s->set = &set_set;
  s->del = &initial_del;
  s->nbytes = &set_nbytes;

  s->f = f;

  os = n0 - b[1];
  l = n - b[0];

  sr = os < l ? 0 : os - l;
  if (os > 0) {
    s->os = os;
    s->sl = 0;
  } else {
    s->sl = -os;
    s->os = 0;
  }
  sr += s->sl;
  s->l = b[0] - sr;
  s->n = b[0];

  return s;
}
