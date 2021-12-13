#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "initial_conditions.h"

void initial_conditions ( int nx, int nt, double x[], double t, double h[], 
  double uh[] )

{
  int i;
  double pi = 3.141592653589793;

  for ( i = 0; i < nx; i++ )
  {
    h[i] = 2.0 + sin ( 2.0 * pi * x[i] );
  }
  for ( i = 0; i < nx; i++ )
  {
    uh[i] = 0.0;
  }
  return;
}
/******************************************************************************/
