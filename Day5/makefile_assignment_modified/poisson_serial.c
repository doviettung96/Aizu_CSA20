#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NX 11
#define NY 11

#include "r8mat_rms.h"
#include "rhs.h"
#include "sweep.h"
#include "timestamp.h"
#include "u_exact.h"
#include "uxxyy_exact.h"

int main(int argc, char *argv[]);

/******************************************************************************/

int main(int argc, char *argv[])

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for POISSON_SERIAL.

  Discussion:

    POISSON_SERIAL is a program for solving the Poisson problem.

    This program runs serially.  Its output is used as a benchmark for
    comparison with similar programs run in a parallel environment.

    The Poisson equation

      - DEL^2 U(X,Y) = F(X,Y)

    is solved on the unit square [0,1] x [0,1] using a grid of NX by
    NX evenly spaced points.  The first and last points in each direction
    are boundary points.

    The boundary conditions and F are set so that the exact solution is

      U(x,y) = sin ( pi * x * y )

    so that

      - DEL^2 U(x,y) = pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )

    The Jacobi iteration is repeatedly applied until convergence is detected.

    For convenience in writing the discretized equations, we assume that NX = NY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 October 2011

  Author:

    John Burkardt
*/
{
  int converged;
  double diff;
  double dx;
  double dy;
  double error;
  double f[NX][NY];
  int i;
  int it;
  int it_max = 1000;
  int j;
  int nx = NX;
  int ny = NY;
  double tolerance = 0.000001;
  double u[NX][NY];
  double u_norm;
  double udiff[NX][NY];
  double uexact[NX][NY];
  double unew[NX][NY];
  double unew_norm;
  double x;
  double y;

  dx = 1.0 / (double)(nx - 1);
  dy = 1.0 / (double)(ny - 1);
  /*
  Print a message.
*/
  timestamp();
  printf("\n");
  printf("POISSON_SERIAL:\n");
  printf("  C version\n");
  printf("  A program for solving the Poisson equation.\n");
  printf("\n");
  printf("  -DEL^2 U = F(X,Y)\n");
  printf("\n");
  printf("  on the rectangle 0 <= X <= 1, 0 <= Y <= 1.\n");
  printf("\n");
  printf("  F(X,Y) = pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )\n");
  printf("\n");
  printf("  The number of interior X grid points is %d\n", nx);
  printf("  The number of interior Y grid points is %d\n", ny);
  printf("  The X grid spacing is %f\n", dx);
  printf("  The Y grid spacing is %f\n", dy);
  /*
  Initialize the data.
*/
  rhs(nx, ny, f);
  /*
  Set the initial solution estimate.
  We are "allowed" to pick up the boundary conditions exactly.
*/
  for (j = 0; j < ny; j++)
  {
    for (i = 0; i < nx; i++)
    {
      if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
      {
        unew[i][j] = f[i][j];
      }
      else
      {
        unew[i][j] = 0.0;
      }
    }
  }
  unew_norm = r8mat_rms(nx, ny, unew);
  /*
  Set up the exact solution.
*/
  for (j = 0; j < ny; j++)
  {
    y = (double)(j) / (double)(ny - 1);
    for (i = 0; i < nx; i++)
    {
      x = (double)(i) / (double)(nx - 1);
      uexact[i][j] = u_exact(x, y);
    }
  }
  u_norm = r8mat_rms(nx, ny, uexact);
  printf("  RMS of exact solution = %g\n", u_norm);
  /*
  Do the iteration.
*/
  converged = 0;

  printf("\n");
  printf("  Step    ||Unew||     ||Unew-U||     ||Unew-Exact||\n");
  printf("\n");

  for (j = 0; j < ny; j++)
  {
    for (i = 0; i < nx; i++)
    {
      udiff[i][j] = unew[i][j] - uexact[i][j];
    }
  }
  error = r8mat_rms(nx, ny, udiff);
  printf("  %4d  %14g                  %14g\n", 0, unew_norm, error);

  for (it = 1; it <= it_max; it++)
  {
    for (j = 0; j < ny; j++)
    {
      for (i = 0; i < nx; i++)
      {
        u[i][j] = unew[i][j];
      }
    }
    /*
  UNEW is derived from U by one Jacobi step.
*/
    sweep(nx, ny, dx, dy, f, u, unew);
    /*
  Check for convergence.
*/
    u_norm = unew_norm;
    unew_norm = r8mat_rms(nx, ny, unew);

    for (j = 0; j < ny; j++)
    {
      for (i = 0; i < nx; i++)
      {
        udiff[i][j] = unew[i][j] - u[i][j];
      }
    }
    diff = r8mat_rms(nx, ny, udiff);

    for (j = 0; j < ny; j++)
    {
      for (i = 0; i < nx; i++)
      {
        udiff[i][j] = unew[i][j] - uexact[i][j];
      }
    }
    error = r8mat_rms(nx, ny, udiff);

    printf("  %4d  %14g  %14g  %14g\n", it, unew_norm, diff, error);

    if (diff <= tolerance)
    {
      converged = 1;
      break;
    }
  }

  if (converged)
  {
    printf("  The iteration has converged.\n");
  }
  else
  {
    printf("  The iteration has NOT converged.\n");
  }
  /*
  Terminate.
*/
  printf("\n");
  printf("POISSON_SERIAL:\n");
  printf("  Normal end of execution.\n");
  printf("\n");
  timestamp();

  return 0;
}

#undef NX
#undef NY
