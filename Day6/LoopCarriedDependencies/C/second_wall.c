/* A gettimeofday routine to give access to the wall
   clock timer on most UNIX-like systems.

   You will need to compile with "-DUNDERSCORE"
   to get this to link with FORTRAN on many systems.
*/

#include <sys/time.h>
/* int gettimeofday(struct timeval *tp, struct timezone *tzp); */
#ifdef UNDERSCORE
#  define mysecond mysecond_
#  define second   second_
#endif

double mysecond ()
{ /* mysecond */
/* 
    struct timeval { long        tv_sec;
                     long        tv_usec;        };
    struct timezone { int   tz_minuteswest;
         int        tz_dsttime;      };
*/

    struct timeval tp;
    struct timezone tzp;
    int i;

    i = gettimeofday(&tp,&tzp);
    return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
} /* mysecond */

double second ()
{ /* second */
    return mysecond();
} /* second */

