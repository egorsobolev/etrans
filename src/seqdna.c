#include <stdio.h>
#include <malloc.h>

#include "types.h"
#include "charge.h"

int readparm(FILE *f, charge_parm_t *p)
{
  int i, fcnt;
  double b;
  for (i = 0; i < 4; ++i) {
    fcnt = fscanf(f, "%lf", &b);
    if (fcnt != 1)
      return -1;
    p->bet[i] = (real_t) b;
  }  
  for (i = 0; i < 16; ++i) {
    fcnt = fscanf(f, "%lf", &b);
    if (fcnt != 1)
      return -1;
    p->lam[i] = (real_t) b;
  }
  return 0;
}

int seqscan(FILE *f, char **seq)
{
  char c, *s;
  int n, i;

  n = 4096;
  s = (char *) malloc(n * sizeof(char));
  *seq = s;
  if (!s) return -1;
  i = 0;
  while (!feof(f)) {
    c = (char) fgetc(f) & 0xDF;
    if (c == 'A' || c == 'T' || c == 'C' || c == 'G') {
      if (i >= n) {
	n += 4096;
	s = (char *) realloc(s, n * sizeof(char));
	*seq = s;
	if (!s) return -1;
      }
      s[i] = c;
      ++i;
    } else if (c != '\n')
      return i;
  }
  return i;
}

int seqdna(charge_parm_t *p, int n, const char *seq, real_t *d, real_t *s)
{
  int k, l, i;

  if (n < 1)
    return 1;

  if (seq[0]  == 'G') k = 0;
  else if (seq[0] == 'A') k = 1;
  else if (seq[0] == 'T') k = 2;
  else if (seq[0] == 'C') k = 3;
  else return 2;
  d[0] = (real_t) p->bet[k];
  l = k << 2;
  for (i = 1; i < n; ++i) {
    if (seq[i]  == 'G') k = 0;
    else if (seq[i] == 'A') k = 1;
    else if (seq[i] == 'T') k = 2;
    else if (seq[i] == 'C') k = 3;
    else return 2;
    d[i] = (real_t) p->bet[k];
    s[i - 1] = (real_t) p->lam[l | k];
    l = k << 2;
  }
  return 0;
}
