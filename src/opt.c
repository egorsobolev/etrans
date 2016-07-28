#include "opt.h"

#include <stdlib.h>
#include <string.h>

#define NOPTION 35

void **argtable_mk(etrans_opt_t *o)
{
  int i = 0;
  void **t;
  t = calloc(NOPTION, sizeof(void *));
  if (!t)
    return NULL;

  t[i++] = o->mdl = arg_str0("C", NULL, NULL, "chain model: OSC, PBD (OSC)");
  t[i++] = o->seqfn = arg_file0("s", NULL, NULL, "chain sequence file");
  /*  t[i++] = o->mg4 = arg_lit0("4", NULL, "4th order of matrix integral");*/
  t[i++] = o->mg = arg_int0("g", NULL, NULL, "order of matrix integral: 1,2,4 (2)");
  t[i++] = o->n0 = arg_int0("0", NULL, NULL, "number of site that is used as the origin ((n-1)/2+1)");
  t[i++] = o->outsh = arg_file0("O", NULL, NULL, "output sheme file");
  t[i++] = o->logfn = arg_file0("l", NULL, NULL, "log file");
  t[i++] = o->init = arg_str0("y", NULL, NULL, "initial state type: 1, U, P, S (1)");
  t[i++] = o->initfn = arg_file0("I", NULL, NULL, "initial state file");
  t[i++] = o->lpfn = arg_file0("w", NULL, NULL, "file of final point set"); 
  t[i++] = o->rst = arg_lit0("c", NULL, "continue from control point");
  t[i++] = o->nsamp = arg_int1("N", NULL, NULL, "number of samples per processor");
  t[i++] = o->tmax = arg_dbl0("t", NULL, NULL, "trajectory time (600)");
  t[i++] = o->h = arg_dbl0("h", NULL, NULL, "time step (0.2)");
  t[i++] = o->no = arg_int0("Q", NULL, NULL, "output step (1)");
  t[i++] = o->nthr = arg_dbl0("H", NULL, NULL, "Norm correction threshold (inf)");
  t[i++] = o->ntm =  arg_dbl0("n", NULL, NULL, "Time interval for norm correction (inf)");
  t[i++] = o->temp = arg_dbl0("T", NULL, NULL, "T - temperature, K (300)");
  t[i++] = o->gamma = arg_dbl0("F", NULL, NULL, "gamma - friction (6.0e-3)");
  t[i++] = o->omegaM2 = arg_dbl0("E", NULL, NULL, "omega0^2 - elastic (1.0e-4)");
  t[i++] = o->sigma = arg_dbl0("S", NULL, NULL, "sigma (1.0)");
  t[i++] = o->omegaB2 = arg_dbl0("X", NULL, NULL, "omegaB^2 - dispersion (0.0)");
  t[i++] = o->rho = arg_dbl0("r", NULL, NULL, "rho (0.5)");
  t[i++] = o->epsilon = arg_dbl0("e", NULL, NULL, "epsilon (.35/4.45)");
  t[i++] = o->lambda = arg_dbl0("U", NULL, NULL, "lambda - electric field");
  t[i++] = o->chi = arg_dbl0("L", NULL, NULL, "chi - coupling constant in quantum equation (0.02)");
  t[i++] = o->mu = arg_dbl0("M", NULL, NULL, "mu - coupling constant in classic equation (0.02)");
  t[i++] = o->prmfn = arg_file0("p", NULL, NULL, "file of quantum parameters");

  t[i++] = o->heattm = arg_dbl0("b", NULL, NULL, "Heating time for PBD chain");
  t[i++] = o->memtm = arg_dbl0("i", NULL, NULL, "Memory time of PDB chain");
  t[i++] = o->gamma_h = arg_dbl0("G", NULL, NULL, "gamma_h - friction in the heating procedure");
  /*
  t[i++] = o->nh = arg_int0("g", NULL, NULL, "number of heating steps (0)");
  t[i++] = o->hh = arg_dbl0("u", NULL, NULL, "free oscilator time step (1e-3 * 2pi/omega0)");
  t[i++] = o->nxt = arg_lit0("z", NULL, "start from free oscilator (else random)");
  */
  t[i++] = o->drp = arg_dbl0("d", NULL, NULL, "time between control points, min (20)");
  t[i++] = o->help  = arg_lit0(NULL, "help", "print this help and exit");
  t[i++] = o->ver  = arg_lit0(NULL, "version", "print version information and exit");
  t[i++] = o->outfn = arg_file1(NULL, NULL, "<output file>", NULL);
  t[i++] = o->end = arg_end(20);

  if (arg_nullcheck(t) != 0) {
    arg_freetable(t, NOPTION);
    free(t);
    return NULL;
  }

  o->tmax->dval[0] = 600.0;
  o->h->dval[0] = 0.2;
  o->no->ival[0] = 1;
  o->temp->dval[0] = 300;
  o->gamma->dval[0] = 6E-3;
  o->omegaM2->dval[0] = 1E-4;
  o->omegaB2->dval[0] = 0.0;
  o->sigma->dval[0] = 1.0;
  o->rho->dval[0] = 0.5;
  o->epsilon->dval[0] = .35 / 4.45;
  o->chi->dval[0] = 0.02;
  o->lambda->dval[0] = 0.0;
  o->drp->dval[0] = 20.0;
  o->mdl->sval[0] = "OSC";
  o->init->sval[0] = "1";
  o->ntm->dval[0] = 1.0 / 0.0;
  o->nthr->dval[0] = 1.0 / 0.0;
  o->mg->ival[0] = 2;

  return t;
}

void argtable_del(void **t)
{
  arg_freetable(t, NOPTION);
  free(t);
}

int options_write(FILE *f, int argc, char **argv)
{
  int i, n;
  if (fwrite(&argc, sizeof(int), 1, f) != 1)
    return -1;
  for (i = 0; i < argc; ++i) {
    n = strlen(argv[i]);
    if (fwrite(&n, sizeof(int), 1, f) != 1)
      return -1;
    if (fwrite(argv[i], sizeof(char), n, f) != n)
      return -1;
  }
  return 0;
}

void options_print(int argc, char **argv)
{
  int i;
  printf(argv[0]);
  for (i = 1; i < argc; ++i)
    printf(" %s", argv[i]);
  printf("\n");
}
