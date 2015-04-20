#include "config.h"
#include "types.h"

/*
  input:
  n - matrix size, int(1)
  br - real part of diagonal, real(n)
  bi - imaginary part of diagonal, real(1)
  cr - real part of superdiagonal, real(n-1)
  ci - imaginary part of superdiagonal, real(n-1)
  dr - real part of right hand vector, real(n)
  di - imaginary part of right hand vector, real(n)
  
  output:
  xr - real part of unknown vector, real(n)
  xi - imaginarty part of unknown vector, real(n)

  descripiton:
  solve system of linear equations using Thomas algorithm

  note:
  dr,di - are destroyed
*/
void thomas(int n, real_t *br, real_t bi, real_t* cr, real_t *ci, real_t *dr, real_t *di, real_t *xr, real_t *xi)
{
  int k, k1, n1;
  real_t b2, mr, mi, yr, yi;

  n1 = n - 1;

  for (k = 0; k < n; ++k) {
    xr[k] = dr[k];
    xi[k] = di[k];
  }
  b2 = br[0]*br[0] + bi * bi;
  dr[0] = br[0] / b2;
  di[0] = bi / b2;
  for (k = 1; k < n; ++k) {
    k1 = k - 1;
    /* m = a(k-1) / b(k-1) */
    mr =  cr[k1]*dr[k1] - ci[k1]*di[k1];
    mi = -ci[k1]*dr[k1] - cr[k1]*di[k1];
    /* d(k) = d(k) - m * d(k-1) */
    xr[k] = dr[k] - mr*xr[k1] + mi*xi[k1];
    xi[k] = di[k] - mi*xr[k1] - mr*xi[k1];
    /* b(k) = b(k) - m * c(k-1); */
    yr = br[k] - mr*cr[k1] + mi*ci[k1];
    yi = bi - mi*cr[k1] - mr*ci[k1];
    b2 = yr*yr + yi*yi;
    dr[k] = yr / b2;
    di[k] = yi / b2;
  }
  /* x(n) = d(n) / b(n); */
  yr = xr[n1]*dr[n1] + xi[n1]*di[n1];
  yi = xi[n1]*dr[n1] - xr[n1]*di[n1];
  xr[n1] = yr;
  xi[n1] = yi;

  for (k = n1-1; k >= 0; k--) {
    k1 = k + 1;
    /* y = d(k) - c(k) * x(k+1); */
    yr = xr[k] - cr[k]*xr[k1] + ci[k]*xi[k1];
    yi = xi[k] - ci[k]*xr[k1] - cr[k]*xi[k1];
    /* x(k) = y / b(k); */
    xr[k] = yr*dr[k] + yi*di[k];
    xi[k] = yi*dr[k] - yr*di[k];
  }
}
