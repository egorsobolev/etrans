#ifndef __ETRANS_FPSET_H
#define __ETRANS_FPSET_H

#include "types.h"

int fpset_create(const char *fn, int n, int n0);
int fpset_open(const char *fn, int n, int n0, int nsamp);
int fpset_write(int n, real_t *x);
void fpset_close();

#endif /*__ETRANS_FPSET_H */
