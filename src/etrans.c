#include "config.h"

#include "osc.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <float.h>

#include <math.h>
#ifdef USE_MKL
# include <mkl.h>
#else
# include <gsl/gsl_cblas.h>

void vdMul(int n, const double *a, const double *b, double *y)
{
   int i;
   for (i = 0; i < n; ++i)
     y[i] = a[i] * b[i];
}
void vdSqrt(int n, const double *a, double *y)
{
   int i;
   for (i = 0; i < n; ++i)
     y[i] = sqrt(a[i]);
}

#endif

#ifdef WIN32
#define unlink _unlink
#else
#include <unistd.h>
#endif

#ifdef MPI
#include <mpi.h>
#define walltime() MPI_Wtime()
#else
#include <time.h>
#define walltime()      ((double) clock() / CLOCKS_PER_SEC)
#endif

#include "etrans.h"
#include "opt.h"

int main(int argc, char **argv)
{
  etrans_opt_t opt;
  void **argtable;
  char *seq, *progname = "etrans", *fntmp;
  int sol_size, m, m1, s, n, nstep, nsamp, exitcode, nerrors, qn, i, heat_nstep;
  size_t fnlen, fcnt;
  real_t h, tmax, heat_h, upr, mu;
  double *rbuf, wtm0, wtm1, elapse, droptime, k1, k2;
  eqdata_t d;
  solution_t rs, ri;
  FILE *f;
  int rank, np;
  charge_eq_t ch;

  argtable = argtable_mk(&opt);
  if (!argtable) {
    fprintf(stderr, "%s: insufficient memory\n", progname);
    exit(-7);

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
    printf("etrans 0.14 20-09-2014\n");
    printf("Authors: Egor Sobolev and Dmitry Tikhonov\n");
    printf("Copyright (C) 2011-2014 Institute of Mathematical Problems of Biology RAS\n");
    exitcode = 0;
    goto exit;
  }
  if (nerrors > 0) {
    arg_print_errors(stderr, opt.end, progname);
    fprintf(stderr, "Try '%s --help' for more information.\n", progname);
    exitcode = -8;
    goto exit;
  }
  if (opt.logfn->count > 0) {
    if (!freopen(opt.logfn->filename[0], "w", stdout)) {
      fprintf(stderr, "%s: can't open log file.\n", progname);
      exitcode = -11;
      goto exit;
    }
  }

#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (!rank)
    printf("Running %s in parallel on %d processor(s).\n", progname, np);
#else
  np = 1;
  rank = 0;
  printf("Running %s in sequential.\n", progname);
#endif
  if (!rank) {
    options_print(argc, argv);
    printf("\n");
  }

  tmax = (real_t) opt.tmax->dval[0];
  h = (real_t) opt.h->dval[0];
  nstep = (int) ceil(tmax / h);
  nsamp = opt.nsamp->ival[0];

  if (opt.prmfn->count > 0) {
    f = fopen(opt.prmfn->filename[0], "r");
    if (!f) {
      fprintf(stderr, "%s: can't open file of quantum parameters.\n", progname);
      exitcode = -13;
      goto exit;
    }
    fcnt = readparm(f);
    fclose(f);
    if (fcnt) {
      fprintf(stderr, "%s: can't read quantum parameters.\n", progname);
      exitcode = -14;
      goto exit;
    }
  }

  f = fopen(opt.seqfn->filename[0], "r");
  if (!f) {
    fprintf(stderr, "%s: can't open sequence file.\n", progname);
    exitcode = -9;
    goto exit;
  }
  n = seqscan(f, &seq);
  fclose(f);
  if (!n) {
    fprintf(stderr, "%s: empty sequence, just quit.\n", progname);
    exitcode = 0;
    goto exit0;
  } else if (n < 0) {
    fprintf(stderr, "%s: insufficient memory\n", progname);
    exitcode = -1;
    goto exit0;
  }

  if (charge_init(&ch, n, (real_t) opt.chi->dval[0])) {
    fprintf(stderr, "%s: insufficient memory\n", progname);
    exitcode = -2;
    goto exit1;
  }
  if (seqdna(n, seq, ch.d, ch.s)) {
    fprintf(stderr, "%s: invalid symbol in sequence\n", progname);
    exitcode = -2;
    goto exit15;
  }
  ch.h_revstep = opt.na->ival[0];
  d.ch = &ch;

  d.n = n;
  d.half = (n - 1) / 2;
  /*  d.sv = (real_t) opt.chi->dval[0];*/
  d.h_revstep = opt.na->ival[0];
  d.x0rnd = (opt.nxt->count == 0);

  d.q_outstep = opt.nq->ival[0] ? opt.nq->ival[0] : nstep;
  d.q_h = h * d.q_outstep;
  d.q_nstep = nstep / d.q_outstep;

  d.osc_outstep = opt.ns->ival[0] ? opt.ns->ival[0] : nstep;
  d.osc_h = h * d.osc_outstep;
  d.osc_nstep = nstep / d.osc_outstep;

  m = n * (d.q_nstep + 1);
  s = n * (d.osc_nstep + 1);
  m1 = m + 3 * (nstep + 1) + s;
  sol_size = 2 * m1;

  fnlen = strlen(opt.outfn->filename[0]);

  d.x = (real_t *) malloc(6 * n * sizeof(real_t) + fnlen + 1);
  if (!d.x) {
    printf("%s: insufficient memory\n", progname);
    exitcode = -1;
    goto exit15;
  }
  d.x0 = d.x + 2 * n;
  fntmp = (char *) (d.x0 + 4 * n);

  strcpy(fntmp, opt.outfn->filename[0]);
  fntmp[fnlen - strlen(opt.outfn->basename[0])] = '$';

  if (opt.initfn->count) {
    f = fopen(opt.initfn->filename[0], "r");
    if (!f) {
      fprintf(stderr, "%s: can't open initial state file.\n", progname);
      exitcode = -17;
      goto exit2;
    }
    exitcode = b0_read(f, n, d.x0);
    fclose(f);
    if (exitcode) {
      fprintf(stderr, "%s: can't read initial state\n", progname);
      exitcode = -17;
      goto exit2;
    }
  } else {
    memset(d.x0, 0, 4 * n * sizeof(real_t));
    d.x0[d.half] = 1.0f;
  }
  /*
    cblas_sscal(n - 1, h, d.s0, 1);
  */
  upr = (real_t) opt.omega0->dval[0];
  mu = (real_t) ((opt.mu->count) ? opt.mu->dval[0] : opt.chi->dval[0]);
  /*
  if (osc_init(&d.osc, n,
	       (real_t) opt.temp->dval[0], (real_t) opt.gamma->dval[0],
	       upr, (real_t) opt.xi->dval[0],
	       (real_t) opt.mu->dval[0])) {
    fprintf(stderr, "%s: insufficient memory\n", progname);
    exitcode = -3;
    goto exit2;
  }
  */
  d.chain = (chain_eq_t *) osc_init(n, (real_t) opt.temp->dval[0], (real_t) opt.gamma->dval[0],
		     upr, (real_t) opt.xi->dval[0], (real_t) opt.mu->dval[0]);
  if (!d.chain) {
    fprintf(stderr, "%s: insufficient memory\n", progname);
    exitcode = -3;
    goto exit2;
  }

  if (solution_allocate(m, nstep, s, &ri, rank ? NULL : &rs, &rbuf)) {
    fprintf(stderr, "%s: insufficient memory\n", progname);
    exitcode = -4;
    goto exit3;
  }

  /* these fields arn't modified */
  if (!rank) {
    rs.nstep = nstep;
    rs.h = h;
    rs.n = n;
    rs.pos = d.q_outstep;
    rs.uos = d.osc_outstep;
  }   
  ri.n = n;
  ri.h = h;
  ri.nstep = nstep;
  ri.pos = d.q_outstep;
  ri.uos = d.osc_outstep;

  /* read file if contunue or reset rs */
  if (opt.rst->count > 0) {
    if (!rank) {
      f = fopen(opt.outfn->filename[0], "rb");
      if (!f) {
	fprintf(stderr, "%s: can't open restart file.\n", progname);
	exitcode = -11;
	goto exit4;
      }
      exitcode = solution_read(f, &ri);
      fclose(f);
      if (exitcode) {
	fprintf(stderr, "%s: can't read solution\n", progname);
	exitcode = -11;
	goto exit4;
      }
      rs.nsamp = ri.nsamp;
      vdMul(sol_size, ri.p, ri.p, rs.p);
      cblas_dscal(m1, (double) rs.nsamp - 1.0, rs.p2, 1);
      cblas_daxpy(m1, (double) rs.nsamp, rs.p, 1, rs.p2, 1);
      cblas_dcopy(m1, ri.p, 1, rs.p, 1);
      cblas_dscal(m1, (double) rs.nsamp, rs.p, 1);
    }
  } else {
    if (!rank) {
      memset(rs.p, 0, sol_size * sizeof(double));
      rs.nsamp = 0;
    }
  }
  if (d.x0rnd) {
    if (!rank)
      printf("Langevin inizialization mode: random\n\n");
  } else {
    if (opt.hh->count > 0)
      heat_h = (real_t) opt.hh->dval[0];
    else
      heat_h = d.chain->period / 1000.0;
    heat_nstep = opt.nh->ival[0];

    d.chain->x0(d.chain, d.x);
    if (!rank) {
      printf("Langevin initialization mode: continue\n\n");
      printf("Langevin initial distibution:\n");
      /*
      if (d.osc.kt > 0.0f) {
	k1 = 1.0f / sqrtf(0.5f * n * d.osc.kt);
	printf(" sqrt(2<u^2>w0^2/kT) = %f\n", k1 * cblas_nrm2(n, d.x, 1) / _sqrt(upr));
	printf(" sqrt(2<v^2>/kT) = %f\n", k1 * cblas_nrm2(n, d.x + n, 1) );
      } else {
	printf(" u[i] = 0, v[i] = 0\n");
      }
      */
      if (heat_nstep > 0) {
	printf("\nheating trajectory:\n");
	printf(" h = %f, t = %f\n", heat_h, heat_h * heat_nstep);
      }
      printf("\n");
      fflush(stdout);
    }
    d.chain->equilibrate(d.chain, heat_h, d.x, heat_nstep, 0);
  }
   

  if (!rank) {
    printf("     #      t,s max|P2-1|    min(h)    max(h)     <u^2>     <v^2>\n");
    fflush(stdout);
  }

  droptime = 60.0 * opt.drp->dval[0];
  i = nsamp;
  qn = (int) (REFERENCE_OPS * droptime * h / (tmax * n) + 0.5);
  if (qn < 1) qn = 1;
  else
    if (qn > i) qn = i;
  wtm0 = walltime();
  while (i) {
    /* reset partial statistics */
    memset(ri.p, 0, sol_size * sizeof(double));
    ri.nsamp = 0;
    /* solve quantum equation qn times */
    if (eq(&d, h, tmax, qn, &ri)) {
      fprintf(stderr, "%s: insufficient memory\n", progname);
      exitcode = -5;
      goto exit4;
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
    MPI_Reduce(ri.p, rbuf, sol_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
    if (!rank) {
      cblas_daxpy(sol_size, 1.0, rbuf, 1, rs.p, 1);
      rs.nsamp += ri.nsamp * np;

      k1 = 1.0 / rs.nsamp; 
      /* calc std */
      if (rs.nsamp > 1) {
	k2 = 1.0 / (rs.nsamp - 1.0);
	vdMul(m1, rs.p, rs.p, ri.p);
	cblas_dscal(m1, -k1 * k2, ri.p, 1);
	cblas_daxpy(m1, k2, rs.p2, 1, ri.p, 1);
	vdSqrt(m1, ri.p, ri.p2);
      } else {
	memset(ri.p2, 0, m1 * sizeof(double));
      }
      /* calc mean */
      cblas_dcopy(m1, rs.p, 1, ri.p, 1);
      cblas_dscal(m1, k1, ri.p, 1);

      ri.nsamp = rs.nsamp;

      f = fopen(fntmp, "wb");
      if (!f) {
	fprintf(stderr, "%s: can't reset output file.\n", progname);
	exitcode = -12;
	goto exit4;
      }
      exitcode = solution_write(f, &ri);
      if (!exitcode)
	exitcode = options_write(f, argc, argv);
      fclose(f);
      if (exitcode) {
	fprintf(stderr, "%s: can't write results.\n", progname);
	exitcode = -12;
	goto exit4;
      }

      unlink(opt.outfn->filename[0]);
      rename(fntmp, opt.outfn->filename[0]);
    }
#ifdef MPI
    MPI_Bcast(&qn, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    wtm1 = walltime();
    elapse += wtm1 - wtm0;
    if (!rank)
      printf("%%%% synctime = %.2g, nsamp = %d, elapse = %.2f\n", wtm1 - wtm0, rs.nsamp, elapse);
    wtm0 = wtm1;
  }
  exitcode = 0;

 exit4:
  free(ri.p);
 exit3:
  d.chain->del(d.chain);
 exit2:
  free(d.x);
 exit15:
  charge_del(&ch);
 exit1:
  free(seq);
 exit0:
#ifdef MPI
  if (exitcode)
    MPI_Abort(MPI_COMM_WORLD, exitcode);
  MPI_Finalize();
#endif
 exit:
  argtable_del(argtable);
  exit(exitcode);
}
