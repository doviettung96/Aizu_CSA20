#ifndef H_SWEEP_H
#define H_SWEEP_H
#define NX 11
#define NY 11
void sweep(int nx, int ny, double dx, double dy, double f[NX][NY],
           double u[NX][NY], double unew[NX][NY]);
#endif
/******************************************************************************/
/*
  Purpose:

   SWEEP carries out one step of the Jacobi iteration.

  Discussion:

    Assuming DX = DY, we can approximate

      - ( d/dx d/dx + d/dy d/dy ) U(X,Y) 

    by

      ( U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1) - 4*U(i,j) ) / dx / dy

    The discretization employed below will not be correct in the general
    case where DX and DY are not equal.  It's only a little more complicated
    to allow DX and DY to be different, but we're not going to worry about 
    that right now.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 October 2011

  Author:

    John Burkardt

  Parameters:

    Input, int NX, NY, the X and Y grid dimensions.

    Input, double DX, DY, the spacing between grid points.

    Input, double F[NX][NY], the right hand side data.

    Input, double U[NX][NY], the previous solution estimate.

    Output, double UNEW[NX][NY], the updated solution estimate.
*/
/******************************************************************************/