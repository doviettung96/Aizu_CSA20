#ifndef H_UXXYY_EXACT_H
#define H_UXXYY_EXACT_H
double uxxyy_exact(double x, double y);
#endif
/******************************************************************************/
/*
  Purpose:

    UXXYY_EXACT evaluates ( d/dx d/dx + d/dy d/dy ) of the exact solution.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 October 2011

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the coordinates of a point.

    Output, double UXXYY_EXACT, the value of 
    ( d/dx d/dx + d/dy d/dy ) of the exact solution at (X,Y).
*/
/******************************************************************************/