#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "r8mat_rms.h"

double r8mat_rms(int nx, int ny, double a[NX][NY])
{
    int i;
    int j;
    double v;

    v = 0.0;

    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            v = v + a[i][j] * a[i][j];
        }
    }
    v = sqrt(v / (double)(nx * ny));

    return v;
}