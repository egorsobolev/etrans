/* intput:
 *   ai(n,m), ar(1), bi(n-1), di(n,m), dr(n,m)
 *   n - matrix dimension
 *   m - number of matrix
 * output:
 *   xi(n,m), xr(n,m)
 */
#include "etrans.h"

#ifdef USE_SINGLE
#define REAL float

void thomas1s(int n, const float *ai, float ar, const float *bi, const float *di, const float *dr, float *xi, float *xr, float *work)
#else
#define REAL double

void thomas1d(int n, const double *ai, double ar, const double *bi, const double *di, const double *dr, double *xi, double *xr, double *work)
#endif
{
	int i, i1, n1;
	REAL *mr, *mi, m2, yr, yi, lr, li, b1;

	n1 = n - 1;
    mr = work;
	mi = mr + n;

	m2 = ai[0] * ai[0] + ar * ar;
	mr[0] = ar / m2;
	mi[0] = ai[0] / m2;

	xr[0] = dr[0];
	xi[0] = di[0];

	for (i = 1; i < n; ++i) {
		i1 = i - 1;

		b1 = bi[i1];

		lr = b1 * mi[i1];
		li = b1 * mr[i1];

		mi[i] = ai[i] - lr * b1;
		mr[i] = ar + li * b1;

		m2 = mi[i] * mi[i] + mr[i] * mr[i];
		mi[i] /= m2;
		mr[i] /= m2;

		xr[i] = dr[i] - lr * xr[i1] + li * xi[i1];
		xi[i] = di[i] - li * xr[i1] - lr * xi[i1];
	}

	yr = xr[n1];
	yi = xi[n1];

	xr[n1] = yr * mr[n1] + yi * mi[n1];
	xi[n1] = yi * mr[n1] - yr * mi[n1];

	for (i = n1 - 1; i >= 0; --i) {
		i1 = i + 1;
		yr = xr[i] + bi[i] * xi[i1];
		yi = xi[i] - bi[i] * xr[i1];
        
		xr[i] = yr * mr[i] + yi * mi[i];
		xi[i] = yi * mr[i] - yr * mi[i];
	}
}
#undef REAL