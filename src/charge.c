#include "charge.h"

#include <stdlib.h>
#include <math.h>

size_t charge_nbytes(int n)
{
  return sizeof(charge_eq_t) + 11 * n * sizeof(real_t) + 5 * n * sizeof(int);
}

charge_eq_t *mk_charge(int n, int n0, charge_parm_t *def)
{
  charge_eq_t *eq;
  eq = (charge_eq_t *) malloc(charge_nbytes(n));
  if (!eq)
    return NULL;
  eq->d = (real_t *) (eq + 1);
  eq->s = eq->d + n;
  eq->s0 = eq->s + n;
  eq->work = eq->s0 + n;
  eq->n = n;
  eq->n0 = n0;
  eq->chi = def->chi;
  eq->lambda = def->lambda;

  return eq;
}

int charge_init(charge_eq_t *eq, etrans_opt_t *o)
{
  int i;

  if (o->chi->count)
    eq->chi = o->chi->dval[0];
  if (o->lambda->count)
    eq->lambda = o->lambda->dval[0];

  for (i = 0; i < eq->n; i++) {
    eq->d[i] += eq->lambda * (i - eq->n0);
  }
  return 0;
}

void charge_del(charge_eq_t *eq)
{
  free(eq);
}

int charge_write(const charge_eq_t *eq, FILE *f)
{
  int n;

  n = 2 * eq->n - 1;

  if (fwrite(&eq->n0, sizeof(int), 1, f) != 1)
    return -1;
  if (fwrite(&eq->chi, sizeof(real_t), 1, f) != 1)
    return -1;
  if (fwrite(&eq->lambda, sizeof(real_t), 1, f) != 1)
    return -1;
  if (fwrite(eq->d, sizeof(real_t), n, f) != n)
    return -1;

  return 0;
}

int charge_read(charge_eq_t *eq, FILE *f)
{
  int n;
  n = 2 * eq->n - 1;

  if (fread(&eq->n0, sizeof(int), 1, f) != 1)
    return -1;
  if (fread(&eq->chi, sizeof(real_t), 1, f) != 1)
    return -1;
  if (fread(&eq->lambda, sizeof(real_t), 1, f) != 1)
    return -1;
  
  if (fread(eq->d, sizeof(real_t), n, f) != n)
    return -1;

  return 0;
}
