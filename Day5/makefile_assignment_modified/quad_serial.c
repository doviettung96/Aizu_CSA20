#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "f.h"
#include "timestamp.h"
#include "cpu_time.h"

int main(int argc, char *argv[]);

/******************************************************************************/

int main(int argc, char *argv[])

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for QUAD_SERIAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 October 2011

  Author:

    John Burkardt
*/
{
  double a = 0.0;
  double b = 10.0;
  double error;
  double exact = 0.49936338107645674464;
  int i;
  int n = 10000000;
  double total;
  double wtime;
  double wtime1;
  double wtime2;
  double x;

  timestamp();
  printf("\n");
  printf("QUAD_SERIAL:\n");
  printf("  C version\n");
  printf("  Estimate the integral of f(x) from A to B.\n");
  printf("  f(x) = 50 / ( pi * ( 2500 * x * x + 1 ) ).\n");
  printf("\n");
  printf("  A        = %f\n", a);
  printf("  B        = %f\n", b);
  printf("  N        = %d\n", n);
  printf("  Exact    = %24.16f\n", exact);

  wtime1 = cpu_time();

  total = 0.0;
  for (i = 0; i < n; i++)
  {
    x = ((n - i - 1) * a + (i)*b) / (n - 1);
    total = total + f(x);
  }

  wtime2 = cpu_time();

  total = (b - a) * total / (double)n;
  error = fabs(total - exact);
  wtime = wtime2 - wtime1;

  printf("\n");
  printf("  Estimate = %24.16f\n", total);
  printf("  Error    = %e\n", error);
  printf("  Time     = %f\n", wtime);
  /*
  Terminate.
*/
  printf("\n");
  printf("QUAD_SERIAL:\n");
  printf("  Normal end of execution.\n");
  printf("\n");
  timestamp();

  return 0;
}
/*******************************************************************************/
