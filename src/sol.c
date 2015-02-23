#include "sol.h"

#include <stdlib.h>
#include <string.h>

int ifunc_get_id(const char *name)
{
  int id;
  id = 0;
  while (id < MAX_FUNC && strcmp(name, int_func[id].name))
    id++;
  if (id >= MAX_FUNC)
    return -1;
  return id;
}

int sol_write_meta(FILE *f, int nfunc, const sol_meta_t *m)
{
  int i;
  unsigned char len;

  if (fwrite(&nfunc, sizeof(int), 1, f) != 1)
    return -1;
  for (i = 0; i < nfunc; i++) {
    len = strlen(int_func[m[i].id].name);
    if (fwrite(&len, sizeof(unsigned char), 1, f) != 1)
      return -1;
    if (fwrite(int_func[m[i].id].name, sizeof(char), len, f) != len)
      return -1;
    if (fwrite(&m[i].flag, sizeof(unsigned char), 1, f) != 1)
      return -1;
    if (m[i].flag & ET_MEAN) {
      if (fwrite(&m[i].mean_step, sizeof(long), 1, f) != 1)
	return -1;
    }
    if (m[i].flag & ET_STD) {
      if (fwrite(&m[i].std_step, sizeof(long), 1, f) != 1)
	return -1;
    }
  }
  return 0;
}

int sol_read_meta(FILE *f, sol_meta_t *m)
{
  int i, nfunc, nb;
  unsigned char len;
  char name[11];
  long n;
  if (fread(&nfunc, sizeof(int), 1, f) != 1)
    return -1;
  for (i = 0; i < nfunc; i++) {
    if (fread(&len, sizeof(unsigned char), 1, f) != 1)
      return -1;
    if (len > 10)
      return -2;
    if (fread(name, sizeof(char), len, f) != len)
      return -1;
    name[len] = '\0';
    m[i].id = ifunc_get_id(name);
    if (m[i].id < 0)
      return -3;
    if (fread(&m[i].flag, sizeof(unsigned char), 1, f) != 1)
      return -1;
    if (m[i].flag & ET_MEAN) {
      if (fread(&m[i].mean_step, sizeof(long), 1, f) != 1)
	return -1;
    } else
      m[i].mean_step = 0;
    if (m[i].flag & ET_STD) {
      if (fread(&m[i].std_step, sizeof(long), 1, f) != 1)
	return -1;
    } else
      m[i].std_step = 0;
  }
  return nfunc;
}

int sol_scan_meta(FILE *f, sol_meta_t *m)
{
  int i, j, id, mean, std, nr;
  char name[11];
  i = 0;
  nr = fscanf(f, "%10s%d%d", name, &mean, &std);
  while (i < MAX_FUNC && !feof(f)) {
    if (nr != 3)
      return -1;
    if (mean < 0 && std < 0)
      return -2;
    id = ifunc_get_id(name);
    if (id < 0)
      return -3;
    m[i].id = id;
    if (mean >= 0) {
      m[i].flag = ET_MEAN;
      m[i].mean_step = mean;
    } else {
      m[i].flag = 0;
      m[i].mean_step = 0;
    }
    if (std >= 0) {
      m[i].flag |= ET_STD;
      m[i].std_step = std;
    } else {
      m[i].std_step = 0;
    }
    i++;
    nr = fscanf(f, "%10s%d%d", name, &mean, &std);
  }
  return i;
}

int sol_default_meta(sol_meta_t *m)
{
  int i;
  char *fn[] = {"x2", "dx2", "d2x2", "Eq", "Eb", "Eu", "Ev", "invb4", "j"};
  for (i = 0; i < 9; i++) {
    m[i].id = ifunc_get_id(fn[i]);
    m[i].flag = ET_MEAN | ET_STD;
    m[i].mean_step = 1;
    m[i].std_step = 1;
  }
  m[i].id = ifunc_get_id("p");
  m[i].flag = ET_MEAN;
  m[i].mean_step = 10;
  m[i].std_step = 0;
  i++;

  return i;
}

void sol_print_meta(FILE *f, int n, const sol_meta_t *m)
{
  int i;
  size_t nb;
  fprintf(f, "  # func     M  mean step S   std step\n");
  for (i = 0; i < n; i++) {
    fprintf(f, "%3d %-8s %c %10d %c %10d\n", i+1, int_func[m[i].id].name,
	    (m[i].flag & ET_MEAN) == ET_MEAN ? '+':'-', m[i].mean_step,
	    (m[i].flag & ET_STD) == ET_STD ? '+':'-', m[i].std_step);
  }
}


sol_t *mk_sol(int nfunc, sol_meta_t *m, int n, long nstep)
{
  int i;
  long numel;
  double *r;
  sol_t *s;
  numel = 0;
  for (i = 0; i < nfunc; i++) {
    if (m[i].flag & ET_MEAN)
      numel += int_func[m[i].id].numel(n, m[i].mean_step ? nstep / m[i].mean_step + 1 : 2);
    if (m[i].flag & ET_STD)
      numel += int_func[m[i].id].numel(n, m[i].std_step ? nstep / m[i].std_step + 1 : 2);
  }
  s = (sol_t *) malloc(sizeof(sol_t) + nfunc * (sizeof(int) + 4*sizeof(long) + 2*sizeof(double *)) + (n + numel) * sizeof(double));
  if (!s)
    return NULL;
  s->imean = (long *) (s + 1);
  s->istd = s->imean + nfunc;
  s->smean = s->istd + nfunc;
  s->sstd = s->smean + nfunc;
  s->nb = (int *) (s->sstd + nfunc);
  s->b = (double *) (s->nb + nfunc);
  s->mean = (double **) (s->b + n);
  s->std = s->mean + nfunc;
  s->r = (double *) (s->std + nfunc);

  s->nsamp = 0;
  s->n = n;
  s->nstep = nstep;
  s->nfunc = nfunc;
  s->m = m;
  s->numel = numel;
  r = s->r;
  for (i = 0; i < nfunc; i++) {
    s->nb[i] = int_func[m[i].id].vlen(n);
    if (m[i].flag & ET_MEAN) {
      s->mean[i] = r;
      if (!m[i].mean_step)
	m[i].mean_step = nstep;
      r += int_func[m[i].id].numel(n, nstep / m[i].mean_step + 1);
      s->smean[i] = 0;
      s->imean[i] = 0;
    } else {
      s->mean[i] = NULL;
    }
    if (m[i].flag & ET_STD) {
      s->std[i] = r;
      if (!m[i].std_step)
	m[i].std_step = nstep;
      r += int_func[m[i].id].numel(n, nstep / m[i].std_step + 1);
      s->sstd[i] = 0;
      s->istd[i] = 0;
    } else {
      s->std[i] = NULL;
    }
  }
  return s;
}

void sol_reset(sol_t *s)
{
  int i;
  for (i = 0; i < s->nfunc; i++) {
    if (s->m[i].flag & ET_MEAN) {
      s->smean[i] = 0;
      s->imean[i] = 0;
    }
    if (s->m[i].flag & ET_STD) {
      s->sstd[i] = 0;
      s->istd[i] = 0;
    }
  }
}

void sol_setzero(sol_t *s)
{
  s->nsamp = 0;
  memset(s->r, 0, s->numel * sizeof(double));
}

void sol_update(const eqdata_t *eq, const real_t *x, const real_t *b2, const real_t *dx, const real_t *d2x, sol_t *s)
{
  int i, j, flag;
  for (i = 0; i < s->nfunc; i++) {
    flag = 1;
    if (s->m[i].flag & ET_MEAN) {
      if (s->smean[i] <= 0) {
	flag = 0;
	int_func[s->m[i].id].f(eq, x, b2, dx, d2x, s->b);
	for (j = 0; j < s->nb[i]; j++)
	  s->mean[i][s->imean[i]++] += s->b[j];
	s->smean[i] = s->m[i].mean_step - 1;
      } else {
	s->smean[i]--;
      }
    }
    if (s->m[i].flag & ET_STD) {
      if (s->sstd[i] <= 0) {
	if (flag)
	  int_func[s->m[i].id].f(eq, x, b2, dx, d2x, s->b);
	for (j = 0; j < s->nb[i]; j++)
	  s->std[i][s->istd[i]++] += s->b[j] * s->b[j];	
	s->sstd[i] = s->m[i].std_step - 1;
      } else {
	s->sstd[i]--;
      }
    }    
  }
}

void sol_normalize(sol_t *s, int nsamp)
{
  long i;
  for (i = 0; i < s->numel; i++)
    s->r[i] /= nsamp;
}
