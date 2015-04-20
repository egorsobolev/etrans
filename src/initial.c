#include "initial.h"
#include "rnd.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

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
  int k, nn = 2 * n;

  for (k = 0; k < nn; ++k) {
    x[k] = s->x0[k];
    x[k+nn] += s->x0[k+nn];
  }

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
  int k, err, m;
#ifdef MPI
  MPI_Status status;
  err = MPI_File_read_shared(s->fh, s->x0, 4 * s->n, MPI_ET_REAL, &status);
  if (err != MPI_SUCCESS)
    return -1;
  err = MPI_Get_count(&status, MPI_ET_REAL, &m);
  if (err != MPI_SUCCESS || m != (4 * s->n))
    return -1;
#else
  m = fread(s->x0, sizeof(real_t), 4 * s->n, s->fh);
  if (m != (4 * s->n))
    return -1;
#endif
  memset(x, 0, 2 * n * sizeof(real_t));

  for (k = 0; k < s->l; ++k) {
    x[s->os+k] = s->x0[s->sl+k];
    x[n+s->os+k] = s->x0[s->n+s->sl+k];
    x[2*n+s->os+k] = s->x0[2*s->n+s->sl+k];
    x[3*n+s->os+k] = s->x0[3*s->n+s->sl+k];
  }

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

void initial_set_del(initial_t *state)
{
  initial_set_t *s = (initial_set_t *) state;
#ifdef MPI
  MPI_File_close(&s->fh);
#else
  fclose(s->fh);
#endif
  free(state);
}

initial_set_t *mk_initial_set(const char *fn, int from, int n, int n0)
{
  initial_set_t *s;
  int b[2], os, sr, l;

#ifdef MPI
  MPI_Status status;
  MPI_File fh;
  if (MPI_File_open(MPI_COMM_WORLD, fn, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh) != MPI_SUCCESS)
    return NULL;
  if (MPI_File_read_at_all(fh, 0, b, 2, MPI_INT, &status) != MPI_SUCCESS)
    return NULL;
  if (MPI_Get_count(&status, MPI_INT, &l) != MPI_SUCCESS || l != 2)
    return NULL;
  if (MPI_File_seek_shared(fh, 2 * sizeof(int) + 4 * b[0] * from * sizeof(real_t), MPI_SEEK_SET) != MPI_SUCCESS)
    return NULL;
#else
  FILE *fh;
  fh = fopen(fn, "r");
  if (!fh)
    return NULL;
  if (fread(&b, sizeof(int), 2, fh) != 2)
    return NULL;
  if (fseek(fh, 4 * b[0] * from, SEEK_CUR))
    return NULL;
#endif

  s = (initial_set_t *) malloc(set_nbytes_impl(b[0]));
  if (!s)
    return NULL;
  s->x0 = (real_t *) (s + 1);

  s->set = &set_set;
  s->del = &initial_set_del;
  s->nbytes = &set_nbytes;

#ifdef MPI
  memcpy(&s->fh, &fh, sizeof(MPI_File));
#else
  s->fh = fh;
#endif

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
