#include <stdlib.h>
#include <stdio.h>
#include "rhs.h"
#include "u_exact.h"
#include "uxxyy_exact.h"
#include "r8mat_rms.h"

void rhs(int nx, int ny, double f[NX][NY])
{
    double fnorm;
    int i;
    int j;
    double x;
    double y;
    /*
  The "boundary" entries of F store the boundary values of the solution.
  The "interior" entries of F store the right hand sides of the Poisson equation.
*/
    for (j = 0; j < ny; j++)
    {
        y = (double)(j) / (double)(ny - 1);
        for (i = 0; i < nx; i++)
        {
            x = (double)(i) / (double)(nx - 1);
            if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
            {
                f[i][j] = u_exact(x, y);
            }
            else
            {
                f[i][j] = -uxxyy_exact(x, y);
            }
        }
    }

    fnorm = r8mat_rms(nx, ny, f);

    printf("  RMS of F = %g\n", fnorm);

    return;
}