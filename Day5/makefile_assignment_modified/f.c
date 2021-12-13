#include "f.h"
double f(double x)
{
    double pi = 3.141592653589793;
    double value;

    value = 50.0 / (pi * (2500.0 * x * x + 1.0));

    return value;
}