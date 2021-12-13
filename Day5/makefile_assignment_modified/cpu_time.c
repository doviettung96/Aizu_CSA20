#include <time.h>
#include "cpu_time.h"
double cpu_time(void)

{
    double value;

    value = (double)clock() / (double)CLOCKS_PER_SEC;

    return value;
}
/******************************************************************************/