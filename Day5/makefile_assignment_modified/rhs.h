#ifndef H_RHS_H
#define H_RHS_H

#define NX 11
#define NY 11
void rhs(int nx, int ny, double f[NX][NY]);
#endif

/******************************************************************************/
/*
  Purpose:

    RHS initializes the right hand side "vector".

  Discussion:

    It is convenient for us to set up RHS as a 2D array.  However, each
    entry of RHS is really the right hand side of a linear system of the
    form

      A * U = F

    In cases where U(I,J) is a boundary value, then the equation is simply

      U(I,J) = F(i,j)

    and F(I,J) holds the boundary data.

    Otherwise, the equation has the form

      (1/DX^2) * ( U(I+1,J)+U(I-1,J)+U(I,J-1)+U(I,J+1)-4*U(I,J) ) = F(I,J)

    where DX is the spacing and F(I,J) is the value at X(I), Y(J) of

      pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 October 2011

  Author:

    John Burkardt

  Parameters:

    Input, int NX, NY, the X and Y grid dimensions.

    Output, double F[NX][NY], the initialized right hand side data.
*/
/******************************************************************************/
