#include <mkl.h>

float clapack_slamch(char cmach)
{
	return slamch(&cmach);
}

void clapack_sstevx(char jobz, char range, int n, float *d, float *e,
			   float vl, float vu, int il, int iu, float abstol,
			   int *m, float *w, float *z, int ldz,
			   float *work, int *iwork, int *ifail, int *info)
{
	sstevx(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, iwork, ifail, info);

}
