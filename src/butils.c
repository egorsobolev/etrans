#include "etrans.h"

#include <memory.h>
#include <malloc.h>
#include <stdio.h>

void diff_b(eqdata_t *d, real_t *b, real_t *u, real_t *db, real_t *d2b)
{
    int n, n1, i;
    real_t *x, *y, *dx, *dy, *d2x, *d2y, *v;
    n = d->n;
    n1 = n - 1;
    v = u + n;
    x = b;
    y = b + n;
    dx = db;
    dy = db + n;
    d2x = d2b;
    d2y = d2b + n;
   
    dx[0] = d->s[0] * y[1] + (d->d[0] + d->sv * u[0]) * y[0];
    dy[0] = -(d->s[0] * x[1] + (d->d[0] + d->sv * u[0]) * x[0]);
    for (i = 1; i < n1; ++i) {
        dx[i] = d->s[i-1] * y[i-1] + d->s[i] * y[i+1] + (d->d[i] + d->sv * u[i]) * y[i];
        dy[i] = -(d->s[i-1] * x[i-1] + d->s[i] * x[i+1] + (d->d[i] + d->sv * u[i]) * x[i]);
    }	
    dx[n1] = d->s[n1-1] * y[n1-1] + (d->d[n1] + d->sv * u[n1]) * y[n1];
    dy[n1] = -(d->s[n1-1] * x[n1-1] + (d->d[n1] + d->sv * u[n1]) * x[n1]);


    d2x[0] = d->s[0] * dy[1] + (d->d[0] + d->sv * u[0]) * dy[0] + d->sv * v[0] * y[0];
    d2y[0] = -(d->s[0] * dx[1] + (d->d[0] + d->sv * u[0]) * dx[0] + d->sv * v[0] * x[0]);
    for (i = 1; i < n1; ++i) {
        d2x[i] = d->s[i-1] * dy[i-1] + d->s[i] * dy[i+1] + (d->d[i] + d->sv * u[i]) * dy[i] + d->sv * v[i] * y[i];
        d2y[i] = -(d->s[i-1] * dx[i-1] + d->s[i] * dx[i+1] + (d->d[i] + d->sv * u[i]) * dx[i] + d->sv * v[i] * x[i]);
    }	
    d2x[n1] = d->s[n1-1] * dy[n1-1] + (d->d[n1] + d->sv * u[n1]) * dy[n1] + d->sv * v[n1] * y[n1];
    d2y[n1] = -(d->s[n1-1] * dx[n1-1] + (d->d[n1] + d->sv * u[n1]) * dx[n1] + d->sv * v[n1] * x[n1]);
}

int b0_read(FILE *f, int n, real_t *x)
{
        int m, c, lm, rm, i, j;
        real_t *a;
        double b;   

        a = x;

        c = fscanf(f, "%d", &m);
        if (c != 1)
                return -1;
        c = abs(m - n);
        lm = c / 2;
        rm = c - lm;

        for (i = 0; i < 4; ++i) {
                if (m > n) {
                        for (j = 0; j < lm; ++j) {
                                c = fscanf(f, "%*f");
                                if (c == EOF) return -1;
                        }
                } else {
                        memset(a, 0, lm * sizeof(real_t));
                        a += lm;
                }
                for (j = 0; j < m; ++j) {
                        c = fscanf(f, "%lf", &b);
                        if (c != 1) return -2;
		        *a = b;
                        ++a;
                }
                if (m > n) {
                        for (j = 0; j < rm; ++j) {
                                c = fscanf(f, "%*f");
                                if (c == EOF) return -1;
                        }
                } else {
                        memset(a, 0, rm * sizeof(real_t));
                        a += rm;
                }
        }
        return 0;
}
