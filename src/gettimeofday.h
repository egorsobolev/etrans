#ifndef __ETRANS_GETTIMEOFDAY_H
#define __ETRANS_GETTIMEOFDAY_H

#ifdef WIN32

#include <time.h>

struct timezone
{
  int  tz_minuteswest;
  int  tz_dsttime;
};

int gettimeofday(struct timeval *tv, struct timezone *tz);

#else
#include <sys/time.h>
#endif

#endif /* __ETRANS_GETTIMEOFDAY_H */
