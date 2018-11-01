#include "config.h"

#include "osc.h"
#include "pbd.h"

#include "initial.h"

#include "sol.h"
#include "fpset.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <float.h>
#include <math.h>
#include <rng.h>

#ifdef WIN32
#include <process.h>
#define unlink _unlink
#define getpid _getpid
typedef int pid_t;
#else
#include <unistd.h>
#endif

#include "gettimeofday.h"

#ifdef MPI
#include <mpi/mpi.h>
#define walltime() MPI_Wtime()
#else
#include <time.h>
#define walltime()      ((double) clock() / CLOCKS_PER_SEC)
#endif

#include "etrans.h"
#include "opt.h"


int eq(eqdata_t *d, real_t h, real_t tmax, int nsamp, int mg, sol_t *res);

int eq_read_head(eqdata_t *d, FILE *f)
{
  if (fread(&d->n, sizeof(int), 1, f) != 1)
    return -1;
  if (fread(&d->nstep, sizeof(long), 1, f) != 1)
    return -1;
  if (fread(&d->nskip, sizeof(int), 1, f) != 1)
    return -1;
  if (fread(&d->h, sizeof(real_t), 1, f) != 1)
    return -1;
  return 0;
}
int etrans_write(const eqdata_t *d, const sol_t *si, int ndone, const double *ss, FILE *f)
{
  if (fwrite(&d->n, sizeof(int), 1, f) != 1)
    return -1;
  if (fwrite(&d->nstep, sizeof(long), 1, f) != 1)
    return -1;
  if (fwrite(&d->nskip, sizeof(int), 1, f) != 1)
    return -1;
  if (fwrite(&d->h, sizeof(real_t), 1, f) != 1)
    return -1;

  if (d->chain->write(d->chain, f))
    return -1;

  if (charge_write(d->ch, f))
    return -1;

  if (fwrite(&ndone, sizeof(int), 1, f) != 1)
    return -1;
  if (sol_write_meta(f, si->nfunc, si->m))
    return -1;
  if (fwrite(si->nb, sizeof(int), si->nfunc, f) != si->nfunc)
    return -1;
  if (fwrite(ss, sizeof(double), si->numel, f) != si->numel)
    return -1;

  return 0;
}

int main(int argc, char **argv)
{
  etrans_opt_t opt;
  void **argtable;
  char *seq, *progname = "etrans", *fntmp;
  int nsamp, exitcode, nerrors, qn, i, heat_nstep, nfunc, ndone, n0, mg;
  long k;
  size_t fnlen, fcnt;
  real_t tmax, heat_h, a;
  double *rbuf, wtm0, wtm1, elapse, droptime;
  eqdata_t d;
  sol_meta_t sol_meta[MAX_FUNC];
  sol_t *si;
  double *ss;
  FILE *f;
  int rank, np;
  int rst, err;
  charge_parm_t charge_parm;
  char mdl[4];
  size_t nb, tmem;
  unsigned int seed;
  struct timeval tm;
  pid_t pid;

  argtable = argtable_mk(&opt);
  if (!argtable) {
    fprintf(stderr, "%s: insufficient memory to allocate argtable.\n", progname);
    exit(-1);
  }
  nerrors = arg_parse(argc, argv, argtable);

  if (opt.help->count > 0) {
    printf("Usage: %s", progname);
    arg_print_syntax(stdout, argtable,"\n");
    arg_print_glossary(stdout, argtable,"  %-20s %s\n");
    exitcode = 0;
    goto exit;
  }
  if (opt.ver->count > 0){
    printf("etrans 1.0 24-02-2015\n");
    printf("Authors: Egor Sobolev and Dmitry Tikhonov\n");
    printf("Copyright (C) 2011-2015 Institute of Mathematical Problems of Biology RAS\n");
    exitcode = 0;
    goto exit;
  }
  if (nerrors > 0) {
    arg_print_errors(stderr, opt.end, progname);
    fprintf(stderr, "Try '%s --help' for more information.\n", progname);
    exitcode = -2;
    goto exit;
  }
  if (opt.logfn->count > 0) {
    if (!freopen(opt.logfn->filename[0], "w", stdout)) {
      fprintf(stderr, "%s: can't open log file.\n", progname);
      exitcode = -3;
      goto exit;
    }
  }

  if (opt.mg->ival[0] == 4) {
    mg = 2;
  } else if (opt.mg->ival[0] == 2) {
    mg = 1;
  } else if (opt.mg->ival[0] == 1) {
    mg = 0;
  } else {
    fprintf(stderr, "%s: invalid value of -g.\n", progname);
    fprintf(stderr, "Try '%s --help' for more information.\n", progname);
    exitcode = -50;
    goto exit;
  }

  pid = getpid();
  gettimeofday(&tm, NULL);
  seed = (int) tm.tv_usec * (int) pid;

#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (!rank)
    printf("Running %s in parallel on %d processor(s).\n", progname, np);
  seed = seed * np + rank;
#else
  np = 1;
  rank = 0;
  printf("Running %s in sequential.\n", progname);
#endif
  if (!rank) {
    options_print(argc, argv);
    printf("\n");
  }
  rng_init(seed);
  /*
  if (rnd_init()) {
    fprintf(stderr, "%s: initialazation of random generator failed.\n", progname);
    exitcode = -4;
    goto exit1;
  }
  */

  d.ntm = (real_t) opt.ntm->dval[0];
  d.nthr = (real_t) opt.nthr->dval[0];

  rst = opt.rst->count;
  if (rst) {
    seq = NULL;

    f = fopen(opt.outfn->filename[0], "rb");
    if (!f) {
      fprintf(stderr, "%s: can't open output/restart file.\n", progname);
      exitcode = -7;
      goto exit2;
    }
    /* read head */
    if (eq_read_head(&d, f)) {
      fprintf(stderr, "%s: can't read header from restart/output file.\n", progname);
      exitcode = -8;
      goto exit2;
    }
    tmax = d.h * d.nstep * d.nskip;
    /* read chain */
    if (fread(mdl, sizeof(char), 3, f) != 3) {
      fprintf(stderr, "%s: can't read chain name from output file.\n", progname);
      exitcode = -8;
      goto exit2;
    }
  } else {
    /* read sequence */
    if (!opt.seqfn->count) {
      fprintf(stderr, "-s <file> or -c must be specified.\n");
      fprintf(stderr, "Try '%s --help' for more information.\n", progname);
      exitcode = -2;
      goto exit2;
    }
    f = fopen(opt.seqfn->filename[0], "r");
    if (!f) {
      fprintf(stderr, "%s: can't open sequence file.\n", progname);
      exitcode = -5;
      goto exit2;
    }
    d.n = seqscan(f, &seq);
    fclose(f);
    if (!d.n) {
      fprintf(stderr, "%s: empty sequence, just quit.\n", progname);
      exitcode = 0;
      goto exit2;
    } else if (d.n < 0) {
      fprintf(stderr, "%s: insufficient memory to allocate sequence.\n", progname);
      exitcode = -6;
      goto exit2;
    }
    /* init head */
    tmax = (real_t) opt.tmax->dval[0];
    d.h = (real_t) opt.h->dval[0];
    d.nskip = opt.no->ival[0];
    d.nstep = (long) ceil(tmax / d.nskip / d.h);
    /* read solution sheme */
    if (opt.outsh->count) {
      f = fopen(opt.outsh->filename[0], "rt");
      if (!f) {
	fprintf(stderr, "%s: can't open file of output sheme.\n", progname);
	exitcode = -9;
	goto exit3;
      }
      nfunc = sol_scan_meta(f, sol_meta);
      if (nfunc < 0) {
	fprintf(stderr, "%s: can't read output sheme from file.\n", progname);
	exitcode = -10;
	goto exit3;
      }
      fclose(f);
    } else {
      nfunc = sol_default_meta(sol_meta);
    }

    strncpy(mdl, opt.mdl->sval[0], 3);
  }
  mdl[3] = '\0';

  nsamp = opt.nsamp->ival[0];

  if (!strcmp(mdl, "OSC")) {
    d.chain = (chain_eq_t *) mk_osc(d.n);
    d.charge_parm = &osc_charge_defs;
  } else if (!strcmp(mdl, "PBD")) {
    d.chain = (chain_eq_t *) mk_pbd(d.n);
    d.charge_parm = &pbd_charge_defs;
  } else
    d.chain = NULL;

  if (!d.chain) {
    fprintf(stderr, "%s: insufficient memory to allocate chain data\n", progname);
    exitcode = -11;
    goto exit2;
  }
  if (rst) {
    /* !!! what do with seq? !!! */
    if (d.chain->read(d.chain, f)) {
      fprintf(stderr, "%s: can't read chain data.\n", progname);
      exitcode = -13;
      goto exit4;
    }
  } else {
    if (d.chain->init(d.chain, &opt)) {
      fprintf(stderr, "%s: initialization of chain failed.\n", progname);
      exitcode = -12;
      goto exit4;
    }
    /* read quantum parameters from file */
    if (opt.prmfn->count > 0) {
      f = fopen(opt.prmfn->filename[0], "rt");
      if (!f) {
	fprintf(stderr, "%s: can't open file of quantum parameters.\n", progname);
	exitcode = -14;
	goto exit4;
      }
      fcnt = readparm(f, &charge_parm);
      fclose(f);
      if (fcnt) {
	fprintf(stderr, "%s: can't read quantum parameters.\n", progname);
	exitcode = -15;
	goto exit4;
      }
      charge_parm.chi = d.charge_parm->chi;
      charge_parm.lambda = d.charge_parm->lambda;
      d.charge_parm = &charge_parm;
    }
  }

  /* make charge equation strucutres */
  if (opt.n0->count) {
    n0 = opt.n0->ival[0];
    if (n0 < 0 || n0 >= d.n) {
      fprintf(stderr, "%s: the origin must be inside the chain.\n", progname);
      exitcode = -16;
      goto exit4;
    }
  } else
    n0 = (d.n - 1) / 2;
  
  d.ch = mk_charge(d.n, n0, d.charge_parm);
  if (!d.ch) {
    fprintf(stderr, "%s: insufficient memory to allocate charge data\n", progname);
    exitcode = -17;
    goto exit4;
  }
  if (rst) {
    /* read charge */
    if (charge_read(d.ch, f)) {
      fprintf(stderr, "%s: can't read charge data\n", progname);
      exitcode = -20;
      goto exit5;
    }
    if (fread(&ndone, sizeof(int), 1, f) != 1) {
      fprintf(stderr, "%s: nsamp read error\n", progname);
      exitcode = -21;
      goto exit5;
    }
    nfunc = sol_read_meta(f, sol_meta);
    if (nfunc <= 0) {
      fprintf(stderr, "%s: output sheme read error\n", progname);
      exitcode = -22;
      goto exit5;
    }
    fseek(f, nfunc * sizeof(int), SEEK_CUR);
  } else {
    /* init charge */
    if (seqdna(d.charge_parm, d.n, seq, d.ch->d, d.ch->s)) {
      fprintf(stderr, "%s: invalid symbol in sequence\n", progname);
      exitcode = -18;
      goto exit5;
    }
    if (charge_init(d.ch, &opt)) {
      fprintf(stderr, "%s: charge init error\n", progname);
      exitcode = -19;
      goto exit5;
    }
    ndone = 0;
  }

  si = mk_sol(nfunc, sol_meta, d.n, d.nstep);
  if (!si) {
    exitcode = -23;
    fprintf(stderr, "%s: insufficient memory to allocate output fuctional\n", progname);
    goto exit5;
  }
  if (!rank) {
    ss = calloc(si->numel, sizeof(double));
    if (!ss) {
      exitcode = -24;
      fprintf(stderr, "%s: insufficient memory to allocate the summation buffer\n", progname);
      goto exit6;
    }
  } else
    ss = NULL;

  if (rst) {
    if (!rank) {
      if (fread(ss, sizeof(double), si->numel, f) != si->numel) {
	fprintf(stderr, "%s: can't read output functionals\n", progname);
	exitcode = -25;
	goto exit7;
      }
    }
    fclose(f);
  } else {
    if (!rank)
      memset(ss, 0, si->numel * sizeof(double));
  }

  if (!strcmp(opt.init->sval[0], "1")) {
    d.s = (initial_t *) mk_initial_1site(d.ch->n0);
  } else if (!strcmp(opt.init->sval[0], "U")) {
    d.s = (initial_t *) mk_initial_uniform();
  } else if (!strcmp(opt.init->sval[0], "P")) {
    if (!opt.initfn->count) {
      fprintf(stderr, "%s: if -y P, you must specify initial state file (-I).\n", progname);
      exitcode = -28;
      goto exit7;
    }
    f = fopen(opt.initfn->filename[0], "r");
    if (!f) {
      fprintf(stderr, "%s: can't open initial state file.\n", progname);
      exitcode = -29;
      goto exit7;
    }
    d.s = (initial_t *) mk_initial_1state(f, d.n, d.ch->n0);
    fclose(f);
  } else if (!strcmp(opt.init->sval[0], "S")) {
    if (!opt.initfn->count) {
      fprintf(stderr, "%s: if -y S, you must specify initial state file (-I).\n", progname);
      exitcode = -28;
      goto exit7;
    }
    d.s = (initial_t *) mk_initial_set(opt.initfn->filename[0], ndone, d.n, d.ch->n0);
  } else {
    fprintf(stderr, "%s: invalid value of -y.\n", progname);
    exitcode = -27;
    goto exit7;
  }
  if (!d.s) {
    exitcode = -26;
    fprintf(stderr, "%s: can't construct initial state structure.\n", progname);
    goto exit8;
  }

  fnlen = strlen(opt.outfn->filename[0]);
  nb = fnlen + 1;
  fntmp = (char *) malloc(nb);
  if (!fntmp) {
    printf("%s: insufficient memory to allocate temporary filename.\n", progname);
    exitcode = -30;
    goto exit9;
  }

  strcpy(fntmp, opt.outfn->filename[0]);
  fntmp[fnlen - strlen(opt.outfn->basename[0])] = '$';

  d.chain->equilibrate(d.chain);

  /*
  goto exit7;
  */

  if (opt.lpfn->count) {
    if (rst) {
      err = fpset_open(opt.lpfn->filename[0], d.ch->n, d.ch->n0, ndone);
    } else {
      err = fpset_create(opt.lpfn->filename[0], d.ch->n, d.ch->n0);
    }
    if (err == -1) {
      fprintf(stderr, "%s: can't open file of final point set.\n", progname);
      exitcode = -31;
      goto exit10;
    } else if (err < 0) {
      if (err == -2) {
	fprintf(stderr, "%s: can't write header of final point set in the file.\n", progname);
	exitcode = -32;
      }
      if (err == -3) {
	fprintf(stderr, "%s: can't read header of final point set from the file.\n", progname);
	exitcode = -33;
      }
      if (err == -4) {
	fprintf(stderr, "%s: number of sites in final point set file does not match the current settings.\n", progname);
	exitcode = -34;
      }
      if (err == -5) {
	fprintf(stderr, "%s: final point set file too short.\n", progname);
	exitcode = -35;
      }
      if (err == -6) {
	fprintf(stderr, "%s: can't set position in final point set file.\n", progname);
	exitcode = -36;
      }
      goto exit11;
    }
  }

  if (!rank) {
    sol_print_meta(stdout, nfunc, sol_meta, d.n, d.nstep);

    nb += eq_nbytes(d.n) + sol_nbytes(d.n, si->numel, si->nfunc) + 
      charge_nbytes(d.n) + d.chain->nbytes(d.n) + d.s->nbytes(d.s);

    printf("Used memory on slave:  %12ld\n", nb);
    printf("Used memory on master: %12ld\n\n", nb + si->numel * sizeof(double));

    printf("Magnus order: %d\n\n", 1 << mg);

    printf("Nsamp: %d/%d\n\n", ndone, ndone + np * nsamp);

    /*               1         2         3         4         5         6         7         8
    /*      12345678901234567890123456789012345678901234567890123456789012345678901234567890 */
    printf("       # ====== progress ======         t max|P2-1| max||A||     P(1)     P(N)\n");
    fflush(stdout);
  }

  droptime = 60.0 * opt.drp->dval[0];
  i = nsamp;
  qn = (int) (REFERENCE_OPS * droptime * d.h / (tmax * d.n) + 0.5);
  if (qn < 1) qn = 1;
  else
    if (qn > i) qn = i;
  wtm0 = walltime();
  while (i) {
    /* reset partial statistics */
    sol_setzero(si);
    /* solve quantum equation qn times */
    if (eq(&d, d.h, tmax, qn, mg, si)) {
      fprintf(stderr, "%s: error in equation solution process\n", progname);
      exitcode = -50;
      goto exit11;
    }
    i -= qn;
    wtm1 = walltime();
    elapse += wtm1 - wtm0;
    qn = (int) (droptime * qn / (wtm1 - wtm0) + 0.5);
    if (qn < 1) qn = 1;
    else
      if (qn > i) qn = i;
    wtm0 = wtm1;
#ifdef MPI
    if (rank) {
      MPI_Reduce(si->r, NULL, si->numel, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&si->nsamp, NULL, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
      MPI_Reduce(MPI_IN_PLACE, si->r, si->numel, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, &si->nsamp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif
    if (!rank) {
      a = (real_t) ndone;
      ndone += si->nsamp;
      for (k = 0; k < si->numel; ++k)
	ss[k] = (a * ss[k] + si->r[k]) / (real_t) ndone;

      f = fopen(fntmp, "wb");
      if (!f) {
	fprintf(stderr, "%s: can't reset output file.\n", progname);
	exitcode = -37;
	goto exit11;
      }
      exitcode = etrans_write(&d, si, ndone, ss, f);
      if (!exitcode)
	exitcode = options_write(f, argc, argv);
      fclose(f);
      if (exitcode) {
	fprintf(stderr, "%s: can't write results.\n", progname);
	exitcode = -38;
	goto exit11;
      }

      unlink(opt.outfn->filename[0]);
      rename(fntmp, opt.outfn->filename[0]);
    }
#ifdef MPI
    MPI_Bcast(&qn, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    wtm1 = walltime();
    elapse += wtm1 - wtm0;
    if (!rank) {
      printf("%%%% synctime = %.2g, nsamp = %d, elapse = %.2f\n", wtm1 - wtm0, ndone, elapse);
      fflush(stdout);
    }
    wtm0 = wtm1;
  }
  exitcode = 0;

 exit11:
  fpset_close();
 exit10:
  free(fntmp);
 exit9:
  d.s->del(d.s);
 exit8:
 exit7:
  if (!rank)
    free(ss);
 exit6:
  free(si);
 exit5:
  charge_del(d.ch);
 exit4:
  d.chain->del(d.chain);
 exit3:
  free(seq);
 exit2:
  /*  rnd_finish();*/
 exit1:
#ifdef MPI
  if (exitcode)
    MPI_Abort(MPI_COMM_WORLD, exitcode);
  MPI_Finalize();
#endif
 exit:
  fflush(stdout);
  argtable_del(argtable);
  exit(exitcode);
}
