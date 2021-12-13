#include <stdio.h>
#include <stdlib.h>

#include "boundary_conditions.h"

/******************************************************************************/

void boundary_conditions ( int nx, int nt, double x[], double t, double h[], 
  double uh[] )
{
  int bc;

  bc = 1;
/*
  Periodic boundary conditions on H and UH.
*/
  if ( bc == 1 )
  {
    h[0]     = h[nx-2];
    h[nx-1]  = h[1];
    uh[0]    = uh[nx-2];
    uh[nx-1] = uh[1];
  }
/*
  Free boundary conditions on H and UH.
*/
  else if ( bc == 2 )
  {
    h[0]     = h[1];
    h[nx-1]  = h[nx-2];
    uh[0]    = uh[1];
    uh[nx-1] = uh[nx-2];
  }
/*
  Reflective boundary conditions on UH, free boundary conditions on H.
*/
  else if ( bc == 3 )
  {
    h[0]     =   h[1];
    h[nx-1]  =   h[nx-2];
    uh[0]    = - uh[1];
    uh[nx-1] = - uh[nx-2];
  }
  return;
}
/******************************************************************************/
