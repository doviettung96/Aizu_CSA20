#include <stdio.h>
#include <sys/param.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/times.h>
#include <time.h>
#include <unistd.h>

#ifdef UNDERSCORE
#  define mysecond mysecond_
#  define second   second_
#endif /* #ifdef UNDERSCORE */

double mysecond ()
{ /* mysecond */
    long sec;
    long clocks_per_second;
    double secx;
    struct tms realbuf;

    times(&realbuf);
    clocks_per_second = sysconf(_SC_CLK_TCK);
    secx = (realbuf.tms_stime + realbuf.tms_utime) / (double)clocks_per_second;
    return ((double) secx);
} /* mysecond */

double second ()
{ /* second */
    return mysecond();
} /* second */
