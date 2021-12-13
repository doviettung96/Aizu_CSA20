#include <math.h>
#include "uxxyy_exact.h"
double uxxyy_exact(double x, double y)
{
    double pi = 3.141592653589793;
    double value;

    value = -pi * pi * (x * x + y * y) * sin(pi * x * y);

    return value;
}