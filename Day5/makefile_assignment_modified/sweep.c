#include <stdlib.h>
#include <stdio.h>

#include "sweep.h"
void sweep(int nx, int ny, double dx, double dy, double f[NX][NY],
           double u[NX][NY], double unew[NX][NY])

{
    int i;
    int j;

    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            if (i == 0 || j == 0 || i == nx - 1 || j == ny - 1)
            {
                unew[i][j] = f[i][j];
            }
            else
            {
                unew[i][j] = 0.25 * (u[i - 1][j] + u[i][j + 1] + u[i][j - 1] + u[i + 1][j] + f[i][j] * dx * dy);
            }
        }
    }
    return;
}