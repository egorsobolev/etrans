#include <stdio.h>
#include <malloc.h>

#include "etrans.h"

static float bet[] = {
  0.00f,	/* G 0 00 */
  6.84f,	/* A 1 01 */
  10.00f,	/* T 2 10 */
  10.00f	/* C 3 11 */
};
static float lam[] = {
  1.276f,	/* GG 00 0000 */
  1.352f,	/* GA 01 0001 */
  2.081f,	/* GT 02 0010 */
  1.671f,	/* GC 03 0011 */
  0.744f,	/* AG 10 0100 */
  0.456f,  /* AA 11 0101 */
  1.595f,	/* AT 12 0110 */
  0.927f,	/* AC 13 0111 */
  1.291f,	/* TG 20 1000 */
  1.307f,	/* TA 21 1001 */
  2.400f,	/* TT 22 1010 */
  1.155f,	/* TC 23 1011 */
  0.638f,	/* CG 30 1100 */
  0.441f,	/* CA 31 1101 */
  1.519f,	/* CT 32 1110 */
  0.623f	/* CC 33 1111 */
};

int readparm(FILE *f)
{
  int i, fcnt;
  for (i = 0; i < 4; ++i) {
    fcnt = fscanf(f, "%f", bet + i);
    if (fcnt != 1)
      return -1;
  }  
  for (i = 0; i < 16; ++i) {
    fcnt = fscanf(f, "%f", lam + i);
    if (fcnt != 1)
      return -1;
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

int seqdna(int n, const char *seq, real_t *d, real_t *s)
{
  int k, l, i;

  if (n < 1)
    return 1;

  if (seq[0]  == 'G') k = 0;
  else if (seq[0] == 'A') k = 1;
  else if (seq[0] == 'T') k = 2;
  else if (seq[0] == 'C') k = 3;
  else return 2;
  d[0] = (real_t) bet[k];
  l = k << 2;
  for (i = 1; i < n; ++i) {
    if (seq[i]  == 'G') k = 0;
    else if (seq[i] == 'A') k = 1;
    else if (seq[i] == 'T') k = 2;
    else if (seq[i] == 'C') k = 3;
    else return 2;
    d[i] = (real_t) bet[k];
    s[i - 1] = (real_t) lam[l | k];
    l = k << 2;
  }
  return 0;
}
