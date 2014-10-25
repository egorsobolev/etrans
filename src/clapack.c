float clapack_slamch(char cmach)
{
  return slamch_(&cmach);
}

void clapack_sstevx(char jobz, char range, int n, float *d, float *e,
		    float vl, float vu, int il, int iu, float abstol,
		    int *m, float *w, float *z, int ldz,
		    float *work, int *iwork, int *ifail, int *info)
{
  sstevx_(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, iwork, ifail, info);
}

double clapack_dlamch(char cmach)
{
  return dlamch_(&cmach);
}

void clapack_dstevx(char jobz, char range, int n, double *d, double *e,
		    double vl, double vu, int il, int iu, double abstol,
		    int *m, double *w, double *z, int ldz,
		    double *work, int *iwork, int *ifail, int *info)
{
  dstevx_(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, iwork, ifail, info);
}
