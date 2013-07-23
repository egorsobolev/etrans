#include "config.h"

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

#include <argtable2.h>

#include "etrans.h"

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

int main(int argc, char **argv)
{
	struct arg_file *opt_seqfn = arg_file1("s", NULL, NULL, "chain sequence file");
	struct arg_file *opt_logfn = arg_file0("l", NULL, NULL, "log file");
	struct arg_file *opt_initfn = arg_file0("I", NULL, NULL, "initial state file");
	struct arg_lit *opt_rst = arg_lit0("c", NULL, "continue from control point");
	struct arg_int *opt_nsamp = arg_int1("N", NULL, NULL, "number of samples per processor");
	struct arg_dbl *opt_tmax = arg_dbl0("t", NULL, NULL, "trajectory time (600)");
	struct arg_dbl *opt_h = arg_dbl0("h", NULL, NULL, "time step (0.2)");
	struct arg_int *opt_na = arg_int0("H", NULL, NULL, "frequency of adaptive step revision (15)");
        struct arg_int *opt_nq = arg_int0("Q", NULL, NULL, "frequency of quantum frames output (1)");
	struct arg_int *opt_ns = arg_int0("f", NULL, NULL, "frequency of classic frames output (5)");
	struct arg_dbl *opt_temp = arg_dbl0("T", NULL, NULL, "T - temperature, K (300)");
	struct arg_dbl *opt_fric = arg_dbl0("F", NULL, NULL, "gamma - friction (6.0e-3)");
	struct arg_dbl *opt_elas = arg_dbl0("E", NULL, NULL, "omega0^2 - elastic (1.0e-4)");
	struct arg_dbl *opt_xi = arg_dbl0("X", NULL, NULL, "xi - dispersion (0.0)");
	struct arg_dbl *opt_lambda = arg_dbl0("L", NULL, NULL, "chi - coupling constant in quantum equation (0.02)");
        struct arg_dbl *opt_mu = arg_dbl0("M", NULL, NULL, "mu - coupling constant in classic equation (0.02)");
	struct arg_file *opt_prm = arg_file0("p", NULL, NULL, "file of quantum parameters");
        struct arg_int *opt_nh = arg_int0("g", NULL, NULL, "number of heating steps (0)");
        struct arg_dbl *opt_hh = arg_dbl0("u", NULL, NULL, "free oscilator time step (1e-3 * 2pi/omega0)");
        struct arg_lit *opt_nxt = arg_lit0("z", NULL, "start from free oscilator (else random)");
	struct arg_dbl *opt_drp = arg_dbl0("d", NULL, NULL, "time between control points, min (20)");
	struct arg_lit *opt_help  = arg_lit0(NULL, "help", "print this help and exit");
	struct arg_lit *opt_ver  = arg_lit0(NULL, "version", "print version information and exit");
	struct arg_file *opt_outfn = arg_file1(NULL, NULL, "<output file>", NULL);
	struct arg_end *end = arg_end(20);
	void *argtable[] = {
		opt_seqfn, opt_logfn, opt_initfn, opt_rst, opt_nsamp, opt_tmax, opt_h,
		opt_na, opt_nq, opt_ns, opt_temp, opt_fric, opt_elas, opt_xi, opt_lambda, 
	        opt_mu, opt_prm, opt_nh, opt_hh, opt_nxt, opt_drp, opt_help, opt_ver,
	        opt_outfn, end
	};
	char *seq, *progname = "etrans", *fntmp;
        int sol_size, m, m1, s, n, nstep, nsamp, nskip, exitcode, nerrors, qn, i, heat_nstep;
	size_t fnlen, fcnt;
	real_t h, tmax, heat_h, upr, mu;
	double *rbuf, wtm0, wtm1, elapse, droptime, k1, k2;
	eqdata_t d;
	solution_t rs, ri;
	FILE *f;
	int rank, np;

	if (arg_nullcheck(argtable) != 0)
	{
		/* NULL entries were detected, some allocations must have failed */
		printf("%s: insufficient memory\n", progname);
		exit(-7);
	}
	/* set any command line default values prior to parsing */
	opt_tmax->dval[0] = 600.0;
	opt_h->dval[0] = 0.2;
	opt_na->ival[0] = 15;
        opt_nq->ival[0] = 1;
	opt_ns->ival[0] = 5;
	opt_temp->dval[0] = 300;
	opt_fric->dval[0] = 6E-3;
	opt_elas->dval[0] = 1E-4;
	opt_xi->dval[0] = 0.0;
	opt_lambda->dval[0] = 0.02;
        /* opt_mu->dval[0] = 0.02; */
	opt_drp->dval[0] = 20.0;
        opt_nh->ival[0] = 0;
	/* Parse the command line as defined by argtable[] */
	nerrors = arg_parse(argc, argv, argtable);

	if (opt_help->count > 0) {
		printf("Usage: %s", progname);
		arg_print_syntax(stdout, argtable,"\n");
		arg_print_glossary(stdout, argtable,"  %-20s %s\n");
		exitcode = 0;
		goto exit;
	}
	/* special case: '--version' takes precedence error reporting */
	if (opt_ver->count > 0){
	        printf("etrans 0.12RC 30-11-2012\n");
		printf("Authors: Egor Sobolev (implementation) and Dmitry Tikhonov (solution)\n");
	        printf("Copyright (C) 2011-2012 Institute of Mathematical Problems of Biology RAS\n");
		exitcode = 0;
		goto exit;
	}
	/* If the parser returned any errors then display them and exit */
	if (nerrors > 0) {
		/* Display the error details contained in the arg_end struct.*/
		arg_print_errors(stdout, end, progname);
		printf("Try '%s --help' for more information.\n", progname);
		exitcode = -8;
		goto exit;
	}
	if (opt_logfn->count > 0) {
		if (!freopen(opt_logfn->filename[0], "w", stdout)) {
			printf("%s: can't open log file.\n", progname);
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

	tmax = (real_t) opt_tmax->dval[0];
	h = (real_t) opt_h->dval[0];
        nskip = opt_nq->ival[0];
	nstep = (int) ceil(tmax / (h * nskip));
	nsamp = opt_nsamp->ival[0];

	if (opt_prm->count > 0) {
		f = fopen(opt_prm->filename[0], "r");
		if (!f) {
			printf("%s: can't open file of quantum parameters.\n", progname);
			exitcode = -13;
			goto exit;
		}
		fcnt = readparm(f);
		fclose(f);
		if (fcnt) {
			printf("%s: can't read quantum parameters.\n", progname);
			exitcode = -14;
			goto exit;
		}
	}

	f = fopen(opt_seqfn->filename[0], "r");
	if (!f) {
		printf("%s: can't open sequence file.\n", progname);
		exitcode = -9;
		goto exit;
	}

	n = seqscan(f, &seq);
	fclose(f);
	if (!n) {
		printf("%s: empty sequence, just quit.\n", progname);
		exitcode = 0;
		goto exit0;
	} else if (n < 0) {
		printf("%s: insufficient memory\n", progname);
		exitcode = -1;
		goto exit0;
	}

	d.n = n;
	d.half = (n - 1) / 2;
	d.sv = (real_t) opt_lambda->dval[0];
        d.h_revstep = opt_na->ival[0];
        d.q_outstep = nskip;
        d.x0rnd = (opt_nxt->count == 0);

	d.osc_outstep = opt_ns->ival[0];
	/* d.osc_nskip = (int) (1.0f / h + .5); */
	d.osc_h = h * d.osc_outstep;
	d.osc_nstep = (int) floor(tmax / d.osc_h);

	m = n * (nstep + 1);
	s = n * (d.osc_nstep + 1);
	m1 = m + 3 * (nstep + 1) + s;
	sol_size = 2 * m1;

	fnlen = strlen(opt_outfn->filename[0]);

	d.x = (real_t *) malloc(8 * n * sizeof(real_t) + fnlen + 1);
	if (!d.x) {
		printf("%s: insufficient memory\n", progname);
		exitcode = -1;
		goto exit1;
	}
	d.d = d.x + 2 * n;	/* n */
	d.s = d.d + n;		/* n - 1 */
	d.x0 = d.s + n;
	fntmp = (char *) (d.x0 + 4 * n);

	strcpy(fntmp, opt_outfn->filename[0]);
	fntmp[fnlen - strlen(opt_outfn->basename[0])] = '$';

	if (seqdna(n, seq, d.d, d.s)) {
		printf("%s: invalid symbol in sequence\n", progname);
		exitcode = -2;
		goto exit2;
	}

	if (opt_initfn->count) {
		f = fopen(opt_initfn->filename[0], "r");
		if (!f) {
			printf("%s: can't open initial state file.\n", progname);
			exitcode = -17;
			goto exit2;
		}
		exitcode = b0_read(f, n, d.x0);
		fclose(f);
		if (exitcode) {
			printf("%s: can't read initial state\n", progname);
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
        upr = (real_t) opt_elas->dval[0];
        mu = (real_t) ((opt_mu->count) ? opt_mu->dval[0] : opt_lambda->dval[0]);
	if (osc_init(&d.osc, n,
		(real_t) opt_temp->dval[0], (real_t) opt_fric->dval[0],
		upr, (real_t) opt_xi->dval[0],
		(real_t) opt_mu->dval[0])) {
			printf("%s: insufficient memory\n", progname);
			exitcode = -3;
			goto exit2;
	}

	if (solution_allocate(m, nstep, s, &ri, rank ? NULL : &rs, &rbuf)) {
		printf("%s: insufficient memory\n", progname);
		exitcode = -4;
		goto exit3;
	}

	/* these fields arn't modified */
	if (!rank) {
		rs.nstep = nstep;
		rs.h = h * nskip;
		rs.n = n;
		rs.o_h = d.osc_h;
		rs.o_nstep = d.osc_nstep;
	}   
	ri.n = n;
	ri.h = h * nskip;
	ri.nstep = nstep;
	ri.o_h = d.osc_h;
	ri.o_nstep = d.osc_nstep;

	/* read file if contunue or reset rs */
	if (opt_rst->count > 0) {
		if (!rank) {
			f = fopen(opt_outfn->filename[0], "rb");
			if (!f) {
				printf("%s: can't open restart file.\n", progname);
				exitcode = -11;
				goto exit4;
			}
			exitcode = solution_read(f, &ri);
			fclose(f);
			if (exitcode) {
				printf("%s: can't read solution\n", progname);
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
	    if (opt_hh->count > 0)
	        heat_h = (real_t) opt_hh->dval[0];
	    else
                heat_h = d.osc.period / 1000.0;
	    heat_nstep = opt_nh->ival[0];

	    osc_x0_rand(&d.osc, d.x);
	    if (!rank) {
	        printf("Langevin initialization mode: continue\n\n");
	        printf("Langevin initial distibution:\n");
	        if (d.osc.kt > 0.0f) {
		    k1 = 1.0f / sqrtf(0.5f * n * d.osc.kt);
		    printf(" sqrt(2<u^2>w0^2/kT) = %f\n", k1 * cblas_nrm2(n, d.x, 1) / _sqrt(upr));
		    printf(" sqrt(2<v^2>/kT) = %f\n", k1 * cblas_nrm2(n, d.x + n, 1) );
		} else {
		    printf(" u[i] = 0, v[i] = 0\n");
		}
	        if (heat_nstep > 0) {
		    printf("\nheating trajectory:\n");
		    printf(" h = %f, t = %f\n", heat_h, heat_h * heat_nstep);
		}
	        printf("\n");
	        fflush(stdout);
	    }
	    osc_equilibrate(&d.osc, heat_h, d.x, heat_nstep);
	}
   

        if (!rank) {
	    printf("     #      t,s max|P2-1|    min(h)    max(h)     <u^2>     <v^2>\n");
	    fflush(stdout);
	}

	droptime = 60.0 * opt_drp->dval[0];
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
			printf("%s: insufficient memory\n", progname);
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
				printf("%s: can't reset output file.\n", progname);
				exitcode = -12;
				goto exit4;
			}
			exitcode = solution_write(f, &ri);
		        if (!exitcode)
		            exitcode = options_write(f, argc, argv);
			fclose(f);
			if (exitcode) {
				printf("%s: can't write results.\n", progname);
				exitcode = -12;
				goto exit4;
			}

			unlink(opt_outfn->filename[0]);
			rename(fntmp, opt_outfn->filename[0]);
		}
#ifdef MPI
		MPI_Bcast(&qn, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
		wtm1 = walltime();
		elapse += wtm1 - wtm0;
		if (!rank)
			printf("%%%% synctime = %.2g, nsamp = %d, elapse = %.2f\n", wtm1 - wtm0, rs.nsamp, elapse);
	}
	exitcode = 0;

exit4:
	free(ri.p);
exit3:
	osc_free(&d.osc);
exit2:
	free(d.x);
exit1:
	free(seq);
exit0:
#ifdef MPI
	if (exitcode)
		MPI_Abort(MPI_COMM_WORLD, exitcode);
	MPI_Finalize();
#endif
exit:
	arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0])); 
	exit(exitcode);
}
