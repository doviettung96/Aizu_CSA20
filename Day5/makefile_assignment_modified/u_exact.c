#include <math.h>
#include "u_exact.h"

double u_exact(double x, double y)
{
    double pi = 3.141592653589793;
    double value;

    value = sin(pi * x * y);

    return value;
}