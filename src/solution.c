#include <stdio.h>
#include <malloc.h>
#include <float.h>
#include <math.h>

#include "etrans.h"

int solution_read(FILE *f, solution_t *rs) 
{
	size_t fcnt;
	int n, sol_size;
	float h;

	fcnt = fread(&n, sizeof(int), 1, f);
	if (fcnt != 1 || rs->n != n) {
		printf("etrans: number of sites in restart is absent or different. stop.\n");
		return -1;
	}
	fcnt = fread(&rs->nsamp, sizeof(int), 1, f);
	if (fcnt != 1) {
		printf("etrans: number of samples in restart is absent. stop.\n");
		return -1;
	}
	fcnt = fread(&n, sizeof(int), 1, f);
	if (fcnt != 1|| rs->nstep != n) {
		printf("etrans: number of steps in restart is absent or different. stop.\n");
		return -1;
	}
	fcnt = fread(&h, sizeof(float), 1, f);
	if (fcnt != 1 || fabsf(rs->h - h) > FLT_EPSILON) {
		printf("etrans: time step for quantum equation in restart is absent or different. stop.\n");
		return -1;
	}
	fcnt = fread(&n, sizeof(int), 1, f);
	if (fcnt != 1|| rs->o_nstep != n) {
		printf("etrans: number of oscilator steps in restart is absent or different. stop.\n");
		return -1;
	}
	fcnt = fread(&h, sizeof(float), 1, f);
	if (fcnt != 1 || fabsf(rs->o_h - h) > FLT_EPSILON) {
		printf("etrans: time step for oscilator in restart is absent or different. stop.\n");
		return -1;
	}

        sol_size = 2 * (rs->n * (rs->nstep + rs->o_nstep + 2) + 3 * rs->nstep + 3);
	fcnt = fread(rs->p, sizeof(double), sol_size, f);
	if (fcnt != sol_size) {
		printf("etrans: cannot read data from restart file.\n");
		return -1;
	}

	return 0;
}


int solution_write(FILE *f, solution_t *rs)
{
	size_t fcnt;
	int sol_size;
   
	fcnt = fwrite(&rs->n, sizeof(int), 1, f);
	if (fcnt != 1) return -1;

	fcnt = fwrite(&rs->nsamp, sizeof(int), 1, f);
	if (fcnt != 1) return -1;

	fcnt = fwrite(&rs->nstep, sizeof(int), 1, f);
	if (fcnt != 1) return -1;

	fcnt = fwrite(&rs->h, sizeof(float), 1, f);
	if (fcnt != 1) return -1;

	fcnt = fwrite(&rs->o_nstep, sizeof(int), 1, f);
	if (fcnt != 1) return -1;

	fcnt = fwrite(&rs->o_h, sizeof(float), 1, f);
	if (fcnt != 1) return -1;

        sol_size = 2 * (rs->n * (rs->nstep + rs->o_nstep + 2) + 3 * rs->nstep + 3);
	fcnt = fwrite(rs->p, sizeof(double), sol_size, f);
	if (fcnt != sol_size) return -1;

	return 0;
}


#ifdef MPI
int solution_allocate(int m, int nstep, int s, solution_t *ri, solution_t *rs, double **rbuf)
{
	int sol_size;

	sol_size = 2 * (m + 3 * (nstep + 1) + s);

	ri->p = (double *) malloc((1 + 2 * (rs != NULL)) * sol_size * sizeof(double));
	if (!ri->p) return -1;

	ri->c = ri->p + m;
        ri->dc = ri->c + nstep + 1;
        ri->d2c = ri->dc + nstep + 1;
	ri->u = ri->d2c + nstep + 1;
	ri->p2 = ri->u + s;
	ri->c2 = ri->p2 + m;
        ri->dc2 = ri->c2 + nstep + 1;
        ri->d2c2 = ri->dc2 + nstep + 1;
	ri->u2 = ri->d2c2 + nstep + 1;

	if (rs) {
		rs->p = ri->p + sol_size;
		rs->c = rs->p + m;
                rs->dc = rs->c + nstep + 1;
                rs->d2c = rs->dc + nstep + 1;
	        rs->u = rs->d2c + nstep + 1;
		rs->p2 = rs->u + s;
		rs->c2 = rs->p2 + m;
                rs->dc2 = rs->c2 + nstep + 1;
                rs->d2c2 = rs->dc2 + nstep + 1;
	        rs->u2 = rs->d2c2 + nstep + 1;
		*rbuf = rs->p + sol_size;
	} else
		*rbuf = NULL;

	return 0;
}
#else
int solution_allocate(int m, int nstep, int s, solution_t *ri, solution_t *rs, double **rbuf)
{
	int sol_size;

	sol_size = 2 * (m + 3 * (nstep + 1) + s);

	ri->p = (double *) malloc(2 * sol_size * sizeof(double));
	if (!ri->p) return -1;

	ri->c = ri->p + m;
        ri->dc = ri->c + nstep + 1;
        ri->d2c = ri->dc + nstep + 1;
	ri->u = ri->d2c + nstep + 1;
	ri->p2 = ri->u + s;
	ri->c2 = ri->p2 + m;
        ri->dc2 = ri->c2 + nstep + 1;
        ri->d2c2 = ri->dc2 + nstep + 1;
	ri->u2 = ri->d2c2 + nstep + 1;

	rs->p = ri->p + sol_size;
	rs->c = rs->p + m;
        rs->dc = rs->c + nstep + 1;
        rs->d2c = rs->dc + nstep + 1;
        rs->u = rs->d2c + nstep + 1;
	rs->p2 = rs->u + s;
	rs->c2 = rs->p2 + m;
        rs->dc2 = rs->c2 + nstep + 1;
        rs->d2c2 = rs->dc2 + nstep + 1;
        rs->u2 = rs->d2c2 + nstep + 1;

	*rbuf = ri->p;

	return 0;
}
#endif /* MPI */
