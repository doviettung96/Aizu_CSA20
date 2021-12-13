#include <stdio.h>
#include <time.h>

#define clocktime clocktime_

double clocktime ()
{ /* clocktime */
    clock_t clock_value;
    return ((double)clock()) / (double)CLOCKS_PER_SEC;
} /* clocktime */

